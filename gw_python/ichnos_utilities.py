import warnings
import geopandas as gpd
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from simplify_lines import simplify_to_linestring_3d
import os

def outline_mesh_to_domain_file(gdf, filename):
    if len(gdf) == 0:
        raise ValueError("GeoDataFrame is empty. Cannot extract domain outline.")

    if len(gdf) > 1:
        warnings.warn(
            "GeoDataFrame has more than one polygon. "
            "Only the first record will be used for domain outline.",
            UserWarning
        )

    geom = gdf.geometry.iloc[0]

    if geom.geom_type != "Polygon":
        raise ValueError("Input geometry must be a single Polygon after dissolve().")

    outer_coords = list(geom.exterior.coords)
    holes = [list(ring.coords) for ring in geom.interiors]

    with open(filename, "w") as f:
        # write outline
        f.write(f"{len(outer_coords)} 1\n")
        for x,y, in outer_coords:
            f.write(f"{x:.3f} {y:.3f}\n")

        # write holes
        for hole_coords in holes:
            f.write(f"{len(hole_coords)}  0\n")
            for x, y in hole_coords:
                f.write(f"{x} {y}\n")


def iwfm_strat_2_ichnos_elev(strat_df, fact=0.3048):
    nlay = int((strat_df.shape[1]-2)/2)

    elev_data = np.zeros((strat_df.shape[0], nlay+1))

    elev_data[:,0] = fact*strat_df["GSE"]
    for ilay in range(0, nlay):
        elev_data[:,ilay+1] = elev_data[:,ilay] - fact*strat_df[[f'A{ilay+1}', f'L{ilay+1}']].to_numpy().sum(axis = 1)

    return elev_data

def domain_top_bot_file(filename,data):
    set1_keys = ['Radius', 'Power', 'XYV']
    set2_keys = ['Npoints', 'Nelem', 'Radius', 'XYV', 'MESH']
    has_set1 = all(key in data for key in set1_keys)
    has_set2 = all(key in data for key in set2_keys)
    # Return True if EITHER set is present
    if not (has_set1 or has_set2):
        raise ValueError("Data does not contain keywords for either CLOUD or MESH2D format.")
    if has_set1 and has_set2:
        raise ValueError("Data contains keywords for BOTH CLOUD and MESH2D formats. Only one is allowed.")

    with open(filename, "w") as f:

        if has_set1:
            f.write("CLOUD\n")
            f.write(f'{data["Radius"]:.3f} {data["Power"]:.3f}\n')
            np.savetxt(f, data['XYV'], fmt='%.3f', delimiter=' ')

        elif has_set2:
            f.write("MESH2D\n")
            f.write(f'{data["Npoints"]:.0f} {data["Nelem"]:.0f} {data["Radius"]:.3f}\n')
            np.savetxt(f, data['XYV'], fmt='%.3f', delimiter=' ')
            np.savetxt(f, data['MESH'], fmt='%.d', delimiter=' ')


def calc_elem_barycenters(nodes, mesh):
    node_lookup = nodes.set_index("ID")[["X", "Y"]]
    barycenters = np.empty((mesh.shape[0], 2))
    is_triangle = (mesh['ND4'] == 0).to_numpy()

    quad_mesh_ids = mesh.index[~is_triangle]

    if len(quad_mesh_ids) > 0:
        quad_node_ids = mesh.loc[quad_mesh_ids, ['ND1', 'ND2', 'ND3', 'ND4']].to_numpy().flatten()
        quad_coords = node_lookup.loc[quad_node_ids].to_numpy()
        quad_barycenters = quad_coords.reshape(-1, 4, 2).mean(axis=1)
        barycenters[~is_triangle] = quad_barycenters

    tri_mesh_ids = mesh.index[is_triangle]
    if len(tri_mesh_ids) > 0:
        tri_node_ids = mesh.loc[tri_mesh_ids, ['ND1', 'ND2', 'ND3']].to_numpy().flatten()

        # Look up ALL necessary coordinates at once
        tri_coords = node_lookup.loc[tri_node_ids].to_numpy()  # (3 * N_tris, 2)

        # Reshape to (N_tris, 3, 2) and calculate the mean across the 3 nodes (axis=1)
        tri_barycenters = tri_coords.reshape(-1, 3, 2).mean(axis=1)
        barycenters[is_triangle] = tri_barycenters

    return barycenters

def find_k_nearest_distances(points, k = 3):
    n_points = points.shape[0]
    if k <= 0:
        raise ValueError("k must be greater than 0.")
    if k >= n_points:
        k = n_points - 1
        print(f"Warning: k was reduced to {k} because the number of points is {n_points}.")

    n_neighbors_to_query = k + 1
    neigh = NearestNeighbors(n_neighbors=n_neighbors_to_query, algorithm='auto', metric='euclidean')
    neigh.fit(points)
    distances, _ = neigh.kneighbors(points)
    k_distances = distances[:, 1:]
    return k_distances

def read_ichnos_traj(filepath: str) -> pd.DataFrame:
    S = {
        "Eid": [],
        "Sid": [],
        "ER": [],  # Exit Result
        "P": [],  # Positions (N x 3 np.array)
        "V": [],  # Velocities (N x 1 np.array)
        "Age": [],  # Ages (N x 1 np.array)
        "Xs": [],
        "Ys": [],
        "Xe": [],
        "Ye": [],
        "Vs": [],
        "Ve": [],
    }

    # Initialize variables for block processing
    positions_list = []
    velocities_list = []
    age_list = []
    current_eid = None
    current_sid = None

    try:
        with open(filepath, 'r') as f:
            for line_num, line in enumerate(f, 1):
                # Clean up the line and skip empty lines or comments
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()

                if parts[0] == '-9':
                    # Expected format: -9 eid sid Exit
                    if len(parts) < 4:
                        raise ValueError(f"Warning: Malformed end line on line {line_num}: {line}. Skipping.")

                    end_eid = parts[1]
                    end_sid = parts[2]
                    exit_value = parts[3]

                    # Check if we accumulated data for this block
                    if positions_list and current_eid == end_eid and current_sid == end_sid:
                        # Create the final numpy arrays
                        p_array = np.array(positions_list, dtype=np.float64)
                        v_array = np.array(velocities_list, dtype=np.float64)
                        age_array = np.array(age_list, dtype=np.float64)

                        # Append data to the DataFrame collection structure
                        S["Eid"].append(current_eid)
                        S["Sid"].append(current_sid)
                        S["ER"].append(exit_value)
                        S["P"].append(p_array)
                        S["V"].append(v_array)
                        S["Age"].append(age_array)
                        S["Xs"].append(p_array[0, 0])
                        S["Ys"].append(p_array[0, 1])
                        S["Xe"].append(p_array[-1, 0])
                        S["Ye"].append(p_array[-1, 1])
                        S["Vs"].append(np.linalg.norm(v_array[0]))
                        ve = np.linalg.norm(v_array[-1])
                        if ve != 0.0:
                            S["Ve"].append(ve)
                        else:
                            if len(v_array) > 1:
                                S["Ve"].append(np.linalg.norm(v_array[-2]))
                            else:
                                S["Ve"].append(0.0)

                    elif positions_list and (current_eid != end_eid or current_sid != end_sid):
                        raise ValueError(
                            f"Error: IDs on end line {end_eid}/{end_sid} do not match accumulated block IDs "
                            f"{current_eid}/{current_sid} on line {line_num}. Block may be corrupted.")

                    # Reset block variables regardless of successful processing
                    positions_list = []
                    velocities_list = []
                    age_list = []
                    current_eid = None
                    current_sid = None

                # --- 2. Process Data Line ---
                else:
                    # Expected format: eid sid x y z vx vy vz age (9 fields)
                    if len(parts) not in [8, 10]:
                        raise ValueError(
                            f"Warning: Malformed data line on line {line_num}: Expected 8 or 10 fields, got {len(parts)}. Skipping.")

                    vel_scalar = False
                    if len(parts) == 8:
                        vel_scalar = True

                    try:
                        eid, sid = parts[1], parts[2]
                        # Convert coordinates/velocities to float
                        data_floats = list(map(float, parts[3:]))
                        x, y, z = data_floats[0:3]
                        age = data_floats[-1]
                        if vel_scalar:
                            v = data_floats[3]
                        else:
                            v = data_floats[3:6]

                        # Initialize block IDs if this is the first data line
                        if current_eid is None:
                            current_eid = eid
                            current_sid = sid
                        elif eid != current_eid or sid != current_sid:
                            # This indicates a missing -9 terminator for the previous block.
                            # For robustness, we skip this line but a real-world parser might
                            # finalize the previous block and start a new one here.
                            print(
                                f"Error: ID change detected on line {line_num} without prior -9 terminator. Expected {current_eid}/{current_sid}, got {eid}/{sid}. Skipping data line.")
                            continue

                        # Accumulate data
                        positions_list.append([x, y, z])
                        velocities_list.append(v)
                        age_list.append(age)
                    except ValueError as e:
                        print(f"Error: Data conversion failed on line {line_num}: {line} -> {e}. Skipping.")
                        continue

    except FileNotFoundError:
        print(f"Error: File not found at path: {filepath}")
        return pd.DataFrame()
    except Exception as e:
        print(f"An unexpected error occurred during file reading: {e}")
        return pd.DataFrame()

    return pd.DataFrame(S)

def streamlines_to_gdf(S,crs,epsilon = 10.0):
    S_simp_lines = S.copy()
    S_simp_lines['geometry'] = S_simp_lines['P'].apply(lambda P: simplify_to_linestring_3d(P, epsilon=20.0))
    S_simp_lines = S_simp_lines[S_simp_lines['geometry'].notnull()].copy()
    S_gdf = gpd.GeoDataFrame(S_simp_lines, geometry='geometry', crs=crs)
    return S_gdf

def read_ichnos_multi_traj(base, niter, nproc, pad=4, show_progress=False):
    frames = []

    for iter in range(niter):
        for i in range(nproc):
            fname = f"{base}_iter_{iter:0{pad}d}_iproc_{i:0{pad}d}.traj"
            if show_progress:
                print(f"Reading file [{fname}]")
            if os.path.exists(fname):
                frames.append(read_ichnos_traj(fname))
            else:
                print(f"Warning: missing file {fname}")

    if len(frames) == 0:
        raise FileNotFoundError("No trajectory files were found.")

    return pd.concat(frames, ignore_index=True)


def map_elem_nodes_to_rows(nodes_df, elem_df):
    # Create mapping NodeID → row position
    id_to_row = (
        nodes_df
        .reset_index()
        .set_index('ID')['index']
    )

    # Columns to map
    nd_cols = ['ND1', 'ND2', 'ND3', 'ND4']

    elem_mapped = elem_df.copy()
    # Replace IDs by row indices
    for col in nd_cols:
        elem_mapped[col] = elem_df[col].map(id_to_row).add(1)

    triangles = elem_df['ND4'] == 0
    elem_mapped.loc[triangles, 'ND4'] = 0

    return elem_mapped

def config_main(xyz_type="CLOUD",method="Euler"):
    data = {
        "Velocity": {
            "XYZType": f"{xyz_type}",
            "ConfigFile": ""
        },

        "Domain": {
            "Outline": "",
            "TopFile": "",
            "BottomFile": "",
            "ProcessorPolys": ""
        },

        "StepConfig": {
            "Method": f"{method}",
            "Direction": 1,
            "StepSize": 50,
            "StepSizeTime": 100000,
            "nSteps": 1,
            "nStepsTime": 0,
            "minExitStepSize": 1
        },

        "StoppingCriteria": {
            "MaxIterationsPerStreamline": 3000,
            "MaxProcessorExchanges": 50,
            "AgeLimit": -1,
            "StuckIter": 10,
            "AttractFile": "",
            "AttractRadius": 30
        },

        "InputOutput": {
            "ParticleFile": "",
            "WellFile": "",
            "OutputFile": "",
            "PrintH5": 0,
            "PrintASCII": 1,
            "ParticlesInParallel": 5000,
            "GatherOneFile": 0
        },
        "Other": {
            "Version": "0.5.07",
            "Nrealizations": 1,
            "nThreads": 1,
            "RunAsThread": 1,
            "OutFreq": 10
        }
    }

    if method == "RK45":
        data["AdaptStep"] = {
            "MaxStepSize": 1000,
            "MinStepSize": 0.1,
            "IncreaseRateChange": 1.5,
            "LimitUpperDecreaseStep": 0.15,
            "Tolerance" : 1
        }
    if method == "PECE":
        data["PECE"] = {
            "Order": 6,
            "Tolerance": 0.2
        }

    return data

def config_vel(xyz_type="CLOUD"):
    data = {
        "Velocity": {
            "Prefix": "",
            "LeadingZeros": 4,
            "Suffix": ".vel",
            "Type": "DETRM",
            "TimeStepFile" : "",
            "TimeInterp": 1,
            "RepeatTime": 0,
            "Multiplier": 1
        },
        "Porosity": {
            "Value": 0.1
        },
        "General": {
            "OwnerThreshold": 0.15,
            "FrequencyStat": 100
        }
    }

    if xyz_type == "MESH2D":
        data["MESH2D"] = {
            "NodeFile": "",
            "MeshFile": "",
            "ElevationFile": "",
            "FaceIdFile": "",
            "Nlayers": 1,
            "INTERP": "ELEMENT"
        }
    if xyz_type == "CLOUD":
        data["MESH2D"] = {
            "Scale": 1,
            "Power": 3,
            "InitDiameter": 3000,
            "InitRatio": 20,
            "GraphPrefix": "",
            "Threshold": 0.1
        }


    return data





