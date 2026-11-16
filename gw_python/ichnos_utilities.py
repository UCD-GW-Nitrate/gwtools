from pathlib import Path
import warnings
import geopandas as gpd
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from simplify_lines import simplify_to_linestring_3d
import os
import configparser
import subprocess
import shlex
import time
from typing import Dict, Any, List, Tuple, Union
from scipy.spatial import cKDTree
from shapely.geometry import Polygon
from math import isfinite

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
        "Ae": []
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
                        S["Ae"].append(age_array[-1])
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
            "Type": "DETRM",
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

def run_ichnos(filename_prefix, main_config_data, vel_config_data, ichnos_exe):
    config_main = configparser.ConfigParser()
    config_main.read_dict(main_config_data)

    main_ini = f"{filename_prefix}_main.ini"
    vel_ini = f"{filename_prefix}_vel.ini"

    main_config_data['Velocity']['ConfigFile'] = vel_ini

    config_vel = configparser.ConfigParser()
    config_vel.read_dict(vel_config_data)

    with open(main_ini, "w") as f:
        config_main.write(f)

    with open(vel_ini, "w") as f:
        config_vel.write(f)

    # shlex.quote ensures correct quoting on Windows/Linux
    cmd = f"{shlex.quote(ichnos_exe)} -c {shlex.quote(main_ini)}"
    print(f"\n▶ Running Ichnos:\n{cmd}\n")
    start_time = time.time()

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        shell=True,
        bufsize=1
    )

    output_log = []  # store for error reporting
    for line in process.stdout:
        print(line, end="")  # live stream
        output_log.append(line)  # save for errors

    process.wait()
    elapsed = time.time() - start_time
    print(f"\n⏱ Ichnos finished in {elapsed:.2f} seconds.")

    cmd = f"{ichnos_exe} -c {filename_prefix}_main.ini"
    with open(os.devnull, 'w') as null:
        subprocess.call(cmd, stdout=null, stderr=null)

    if process.returncode != 0:
        full_output = "".join(output_log)
        raise RuntimeError(
            f"Ichnos execution failed with return code {process.returncode}.\n"
            f"Output:\n{full_output}"
        )

    return True


def read_ichnos_elev(filename: str) -> Dict[str, Any]:
    Elev: Dict[str, Any] = {}
    try:
        with open(filename, 'r') as f:
            # 1. Create a line iterator for clean, sequential reading
            lines_iterator = iter(f.readlines())

            # --- Read the interpolation type (Elev.Type) ---
            try:
                # Read the first line
                current_line = next(lines_iterator).strip()
                while not current_line and lines_iterator:
                    current_line = next(lines_iterator).strip()
                if not current_line:
                    raise ValueError(f"File {filename} is empty or contains no data lines.")

                Elev['Type'] = current_line.split()[0].upper()  # Read the first word and capitalize
            except StopIteration:
                raise ValueError("Could not read interpolation type from file.")

            # --- Handle different types ---
            if Elev['Type'] == 'CLOUD':
                # --- Read parameters R and P ---
                try:
                    current_line = next(lines_iterator).strip()
                    while not current_line and lines_iterator:
                        current_line = next(lines_iterator).strip()

                    params = [float(x) for x in current_line.split()]
                    if len(params) < 2:
                        raise ValueError("CLOUD type requires at least 2 parameters (R, P).")
                    Elev['R'] = params[0]
                    Elev['P'] = params[1]
                except StopIteration:
                    # R, P data was expected but file ended
                    raise ValueError("File ended unexpectedly while reading CLOUD parameters (R, P).")
                except ValueError as e:
                    # Non-float data or split error
                    raise ValueError(f"Error parsing CLOUD parameters (R, P): {e}")

                # --- Read Data points (Elev.Data) ---
                pp: List[List[float]] = []
                while True:
                    try:
                        current_line = next(lines_iterator).strip()
                        # Stop if line is empty (matches MATLAB's `isempty(lines{idx,1})`)
                        if not current_line:
                            break

                        # Read all floats in the line
                        data_row = [float(x) for x in current_line.split()]
                        pp.append(data_row)

                    except StopIteration:
                        # File ended (matches MATLAB's `if idx > length(lines)`)
                        break
                    except ValueError as e:
                        # Skip or report invalid data lines if necessary
                        warnings.warn(f"Skipping line with invalid float data: '{current_line.strip()}'. Error: {e}")

                # Convert list of lists to a NumPy array for efficiency and matrix-like structure
                Elev['Data'] = np.array(pp)
            elif Elev['Type'] == 'MESH2D':
                warnings.warn('Reading MESH2D is not implemented yet in this function.')
            else:
                warnings.warn(f"Unknown elevation type: {Elev['Type']}. Data reading skipped.")

    except FileNotFoundError:
        print(f"Error: File not found at path: {filename}")
        raise

    return Elev


def _prepare_ichnos_mesh(xy, mesh):
    """Validate the mesh and normalize connectivity to four CCW columns."""
    xy = np.asarray(xy, dtype=float)
    mesh = np.asarray(mesh)

    if xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all():
        raise ValueError("xy must contain finite coordinates with shape (n_nodes, 2)")

    if mesh.ndim != 2 or mesh.shape[1] not in (3, 4) or len(mesh) == 0:
        raise ValueError("mesh must have shape (n_elements, 3) or (n_elements, 4)")

    if not np.issubdtype(mesh.dtype, np.integer):
        raise ValueError("mesh must contain integer node indices")

    if np.any(mesh[:, :3] < 0) or np.any(mesh >= len(xy)):
        raise ValueError("mesh contains an out-of-range node index")
    if mesh.shape[1] == 4 and np.any(mesh[:, 3] < -1):
        raise ValueError("use -1 in the fourth column to mark triangles")

    conn = np.full((len(mesh), 4), -1, dtype=np.int64)
    conn[:, :mesh.shape[1]] = mesh

    for n_vertices, mask in ((3, conn[:, 3] == -1), (4, conn[:, 3] >= 0)):
        ids = conn[mask, :n_vertices]
        if np.any(np.diff(np.sort(ids, axis=1), axis=1) == 0):
            raise ValueError("an element repeats a node index")

        vertices = xy[ids]
        # Shift before computing signed area to reduce cancellation.
        local = vertices - vertices[:, :1, :]
        x = local[:, :, 0]
        y = local[:, :, 1]
        area2 = np.sum(
            x * np.roll(y, -1, axis=1) - y * np.roll(x, -1, axis=1),
            axis=1,
        )
        if np.any(area2 == 0):
            raise ValueError("mesh contains a zero-area element")

        clockwise = area2 < 0
        # Keep the first vertex and reverse the remaining vertices.
        ids[clockwise, 1:] = ids[clockwise, 1:][:, ::-1]
        conn[mask, :n_vertices] = ids

    return xy, conn


def _ichnos_radius_table(xy, conn, n_neighbors, cover_elements):
    if (
        isinstance(n_neighbors, (bool, np.bool_))
        or not isinstance(n_neighbors, (int, np.integer))
        or not 1 <= n_neighbors <= len(conn)
    ):
        raise ValueError("n_neighbors must be an integer from 1 to n_elements")

    n_points = np.arange(1, n_neighbors + 1)
    is_quad = conn[:, 3] >= 0

    vertex_ids = conn.copy()
    vertex_ids[~is_quad, 3] = conn[~is_quad, 0]
    vertices = xy[vertex_ids]

    barycenters = vertices[:, :3].mean(axis=1)
    barycenters[is_quad] = vertices[is_quad].mean(axis=1)

    tree = cKDTree(barycenters)
    distances, neighbor_ids = tree.query(barycenters, k=n_points.tolist())
    required = distances.max(axis=0)

    if cover_elements:
        required = np.zeros(n_neighbors, dtype=float)
        for element_vertices, ids in zip(vertices, neighbor_ids):
            offsets = element_vertices[:, None, :] - barycenters[ids]
            farthest = np.hypot(offsets[..., 0], offsets[..., 1]).max(axis=0)
            required = np.maximum(required, np.maximum.accumulate(farthest))

    return pd.DataFrame({"n_points": n_points, "_dist": required})


def ichnos_mesh2d_search_radius(
    xy, mesh, n_neighbors=3, *, cover_elements=False
):
    """Return search distances for finding 1 through n_neighbors barycenters."""
    xy, conn = _prepare_ichnos_mesh(xy, mesh)
    return _ichnos_radius_table(xy, conn, n_neighbors, cover_elements)


def write_ichnos_elev_MESH2D(
    filename, xy, mesh, elev, radius=None, *, n_neighbors=3
):
    """Write one elevation surface or combined top and bottom surfaces."""
    xy, conn = _prepare_ichnos_mesh(xy, mesh)

    elev = np.asarray(elev, dtype=float)
    if elev.ndim == 1:
        elev = elev[:, None]

    if elev.shape not in ((len(xy), 1), (len(xy), 2)):
        raise ValueError("elev must have one or two columns and one row per node")
    if not np.isfinite(elev).all():
        raise ValueError("elev must contain finite values")

    if radius is None:
        radius_table = _ichnos_radius_table(
            xy, conn, n_neighbors, cover_elements=True
        )
        radius = float(np.nextafter(radius_table["_dist"].iloc[-1], np.inf))
    else:
        radius = float(radius)
        if not np.isfinite(radius) or radius <= 0:
            raise ValueError("radius must be finite and positive")

    # Zero-based node indices become one-based. Triangle padding becomes 0.
    with Path(filename).expanduser().open("w", encoding="utf-8") as stream:
        stream.write("MESH2D\n")
        stream.write(f"{len(xy)} {len(conn)} {radius:.17g}\n")
        np.savetxt(stream, np.column_stack((xy, elev)), fmt="%.17g")
        np.savetxt(stream, conn + 1, fmt="%d")

    return radius


def write_ichnos_wells(filename, Eid, XY, T, B, RT=None, randomize=False,*, header=None):
    """Write an ICHNOS particle-release well file."""
    settings = {"npart": 100, "nlay": 25, "rad": 10.5}
    if header is not None:
        unknown = header.keys() - settings.keys()
        if unknown:
            raise ValueError(f"Unknown header keys: {sorted(unknown)}")
        settings.update(header)

    ids = np.asarray(Eid)
    if ids.ndim != 1:
        raise ValueError("Eid must be a one-dimensional array")

    ids = pd.to_numeric(pd.Series(ids), errors="raise")
    if ids.isna().any() or not np.isfinite(ids).all() or (ids % 1 != 0).any():
        raise ValueError("Eid must contain finite integer IDs")
    if ids.duplicated().any():
        raise ValueError("Eid must be unique for each well")

    n_wells = len(ids)
    xy = np.asarray(XY, dtype=float)
    top = np.asarray(T, dtype=float)
    bottom = np.asarray(B, dtype=float)

    if xy.shape != (n_wells, 2):
        raise ValueError("XY must have shape (n_wells, 2)")
    if top.shape != (n_wells,) or bottom.shape != (n_wells,):
        raise ValueError("T and B must each contain one elevation per well")

    release_time = np.asarray(0.0 if RT is None else RT, dtype=float)
    if release_time.ndim == 0:
        release_time = np.full(n_wells, release_time.item())
    if release_time.shape != (n_wells,):
        raise ValueError("RT must be a scalar or contain one value per well")

    values = np.column_stack((xy, top, bottom, release_time))
    if not np.isfinite(values).all():
        raise ValueError("XY, T, B, and RT must contain finite values")

    order = np.arange(n_wells)
    if randomize:
        np.random.default_rng().shuffle(order)

    with Path(filename).expanduser().open("w", encoding="utf-8") as f:
        f.write("100 25 10.500\n")
        for i in order:
            x, y, t, b, rt = values[i]
            f.write(f"{int(ids.iat[i])} {x:.3f} {y:.3f} {t:.3f} {b:.3f} {rt:.3f}\n")


def read_ichnos_outline(filename, crs=None):
    """Read ICHNOS outline polygons into a GeoDataFrame."""
    path = Path(filename).expanduser()
    records = []

    with path.open(encoding="utf-8-sig") as stream:
        # Skip blank lines between or within polygon blocks.
        lines = (
            (line_no, line.split())
            for line_no, line in enumerate(stream, start=1)
            if line.strip()
        )

        for header_line, fields in lines:
            try:
                n_vertices, inout = map(int, fields)
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{header_line}: expected 'n_vertices INOUT'"
                ) from exc

            if n_vertices < 3 or inout not in (0, 1):
                raise ValueError(
                    f"{path}:{header_line}: expected at least three vertices "
                    "and INOUT equal to 0 or 1"
                )

            xy = np.empty((n_vertices, 2), dtype=float)
            for i in range(n_vertices):
                row = next(lines, None)
                if row is None:
                    raise ValueError(
                        f"{path}:{header_line}: expected {n_vertices} vertices, "
                        f"found {i} before end of file"
                    )

                line_no, fields = row
                try:
                    x, y = map(float, fields)
                except ValueError as exc:
                    raise ValueError(
                        f"{path}:{line_no}: expected two coordinates, X Y"
                    ) from exc
                if not np.isfinite([x, y]).all():
                    raise ValueError(f"{path}:{line_no}: coordinates must be finite")
                xy[i] = x, y

            polygon = Polygon(xy)
            if not polygon.is_valid or polygon.area == 0:
                raise ValueError(
                    f"{path}:{header_line}: invalid or zero-area polygon"
                )

            records.append(
                {"INOUT": inout, "n_vertices": n_vertices, "geometry": polygon}
            )

    return gpd.GeoDataFrame(
        records,
        columns=["INOUT", "n_vertices", "geometry"],
        geometry="geometry",
        crs=crs,
    )

def read_ichnos_processor_polys(filename, crs=None):
    """Read an ICHNOS ProcessorPolys or ExpandedPolys file."""
    path = Path(filename).expanduser()
    records = []
    seen_ids = set()

    with path.open(encoding="utf-8-sig") as stream:
        lines = (
            (line_no, line.split())
            for line_no, line in enumerate(stream, start=1)
            if line.strip()
        )

        first = next(lines, None)
        if first is None:
            raise ValueError(f"{path}: missing Nproc header")

        line_no, fields = first
        try:
            nproc, = map(int, fields)
        except ValueError as exc:
            raise ValueError(f"{path}:{line_no}: expected one integer for Nproc") from exc

        if nproc < 1:
            raise ValueError(f"{path}:{line_no}: Nproc must be positive")

        for _ in range(nproc):
            header = next(lines, None)
            if header is None:
                raise ValueError(
                    f"{path}: expected {nproc} polygons, found {len(records)}"
                )

            header_line, fields = header
            try:
                proc_id, n_vertices = map(int, fields)
            except ValueError as exc:
                raise ValueError(f"{path}:{header_line}: expected 'ID n_vertices'") from exc

            if not 0 <= proc_id < nproc or proc_id in seen_ids:
                raise ValueError(
                    f"{path}:{header_line}: duplicate or out-of-range ID {proc_id}"
                )
            if n_vertices < 3:
                raise ValueError(f"{path}:{header_line}: expected at least three vertices")

            xy = []
            for _ in range(n_vertices):
                row = next(lines, None)
                if row is None:
                    raise ValueError(
                        f"{path}:{header_line}: processor {proc_id} declares "
                        f"{n_vertices} vertices, found {len(xy)} before end of file"
                    )

                line_no, fields = row
                try:
                    x, y = map(float, fields)
                except ValueError as exc:
                    raise ValueError(f"{path}:{line_no}: expected 'X Y'") from exc

                if not (isfinite(x) and isfinite(y)):
                    raise ValueError(f"{path}:{line_no}: coordinates must be finite")
                xy.append((x, y))

            polygon = Polygon(xy)
            if not polygon.is_valid or polygon.area == 0:
                raise ValueError(f"{path}:{header_line}: invalid polygon for ID {proc_id}")

            seen_ids.add(proc_id)
            records.append(
                {"ID": proc_id, "n_vertices": n_vertices, "geometry": polygon}
            )

        extra = next(lines, None)
        if extra is not None:
            raise ValueError(f"{path}:{extra[0]}: extra data after {nproc} polygons")

    return gpd.GeoDataFrame(records, geometry="geometry", crs=crs)







