import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon

def read_urf_calc(prefix_filename, nproc, n_strml, por):
    col_names = ['SrcID', 'SrcInd', 'Eid', 'Sid', 'ER',
                 'p_cdsX', 'p_cdsY', 'p_cdsZ', 'v_cds',
                 'p_lndX', 'p_lndY', 'Len']

    urf_data = pd.DataFrame(np.nan, index=range(n_strml), columns=col_names)
    urf_msa = np.full((n_strml, 3 * len(por)), np.nan)

    ida = 0
    for ii in range(nproc):
        fname = prefix_filename + f"{ii}.dat"
        #print("update2")
        df = pd.read_table(fname, sep=',', index_col=False)
        df.columns = df.columns.str.strip()

        nrows = df.shape[0]
        idb = ida + nrows
        urf_data.iloc[ida:idb, urf_data.columns.get_loc('SrcID')] = int(ii)
        urf_data.iloc[ida:idb, urf_data.columns.get_loc('SrcInd')] = range(ida, idb)
        for col in ['Eid', 'Sid', 'ER', 'v_cds',
                    'p_cdsX', 'p_cdsY', 'p_cdsZ',
                    'p_lndX', 'p_lndY', 'Len']:
            #print(col)
            if col in df.columns:
                urf_data.iloc[ida:idb, urf_data.columns.get_loc(col)] = df[col].values
            else:
                raise KeyError(f"Column '{col}' missing in file {fname}")



        idx = 0
        for k in por:
            for suffix in ['mean', 'std', 'Age']:
                colname = f"{suffix}{k}"
                if colname not in df:
                    raise KeyError(f"Column '{colname}' missing in {fname}")
                urf_msa[ida:idb, idx] = df[colname].to_numpy()
                idx += 1

        ida = idb

    valid_mask = ~urf_data['SrcID'].isna()
    urf_data = urf_data.loc[valid_mask].reset_index(drop=True)
    urf_msa = urf_msa[valid_mask.to_numpy()]

    return urf_data, urf_msa

def calculate_well_travel_time(wells, urf_data, urfs_msa, age_col=2):
    results = []

    # Pre-extract travel times for speed
    ages = urfs_msa[:, age_col]

    for eid in wells['Eid']:
        mask = urf_data['Eid'] == eid
        if not mask.any():
            # If no URF records exist for this well
            results.append({'Eid': eid, 'Age': np.nan})
            continue

        # Extract weights and ages
        weights = urf_data.loc[mask, 'v_cds'].values
        ages_subset = ages[mask.values]

        # Compute weighted average (weights are not normalized)
        total_weight = weights.sum()
        if total_weight == 0:
            avg_age = np.nan
        else:
            avg_age = np.dot(weights, ages_subset) / total_weight

        results.append({'Eid': eid, 'Age': avg_age})

    return pd.DataFrame(results)

def write_scattered(filename, data, tri, hor_interp, ver_interp = None, dim = '3D'):
    n_points, n_cols = data.shape
    n_tri = tri.shape[0]
    n_props = n_cols - 2

    with open(filename, 'w') as f:
        f.write("SCATTERED\n")
        f.write(f"{dim}\n")
        if dim == '3D':
            f.write(f"{hor_interp} {ver_interp}\n")
        else:
            f.write(f"{hor_interp}\n")

        f.write(f"{n_points} {n_props} {n_tri}\n")

        # Write data block
        np.savetxt(f, data, fmt="%.3f")

        # Write triangulation block (convert to 1-based indices for output)
        np.savetxt(f, tri, fmt="%d")

def mesh_to_gdf(mesh, crs):
    nodes = mesh['nodes'][:,:2]
    elements = mesh['elements']

    polygons = []
    for elem in elements:
        try:
            coords = [tuple(nodes[i]) for i in elem]
            poly = Polygon(coords)
            if poly.is_valid:
                polygons.append(poly)
            else:
                polygons.append(poly.buffer(0))
        except Exception as e:
            print(f"Error processing element {elem}: {e}")


    gdf = gpd.GeoDataFrame({'elem_id': range(len(polygons))},
                           geometry=polygons, crs=crs)

    return gdf

def read_streams(filename, crs=3310):
    """
    Reads a custom file with N features, each having nv points and attributes R and W,
    and returns a GeoDataFrame with LineString geometries and R, W fields.
    """
    geometries = []
    R_list = []
    W_list = []

    with open(filename, "r") as f:
        # First line: number of features
        N = int(f.readline().strip())

        for _ in range(N):
            # Read feature header: nv, R, W
            header = f.readline().strip().split()
            nv = int(header[0])
            R = float(header[1])
            if nv == 2:
                W = float(header[2])
            else:
                W = np.nan

            # Read nv lines of coordinates
            coords = [tuple(map(float, f.readline().strip().split())) for _ in range(nv)]
            # coords = []
            # for _ in range(nv):
            #     x, y = map(float, f.readline().strip().split())
            #     coords.append((x, y))
            if nv == 2:
                # two points -> rectangular polygon
                (x1, y1), (x2, y2) = coords
                # direction vector
                dx, dy = x2 - x1, y2 - y1
                length = np.hypot(dx, dy)
                # normalized perpendicular vector
                nx, ny = -dy / length, dx / length
                # offsets (W on each side)
                p1_left = (x1 + W * nx, y1 + W * ny)
                p1_right = (x1 - W * nx, y1 - W * ny)
                p2_left = (x2 + W * nx, y2 + W * ny)
                p2_right = (x2 - W * nx, y2 - W * ny)
                # order the corners to make a proper polygon
                geom = Polygon([p1_left, p2_left, p2_right, p1_right])
            else:
                # more than two points -> polygon from coordinates
                geom = Polygon(coords)

            geometries.append(geom)
            R_list.append(R)
            W_list.append(W)

    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame({"R": R_list, "W": W_list}, geometry=geometries)
    gdf.set_crs(crs, inplace=True)

    return gdf

def read_mesh(filename):
    with open(filename, 'r') as f:
        # --- Read Nnodes and Nelem ---
        first_line = f.readline().strip().split()
        if len(first_line) < 2:
            raise ValueError("First line must contain Nnodes and Nelem.")
        Nnodes, Nelem = map(int, first_line[:2])

        # --- Read node coordinates ---
        nodes = np.loadtxt(f, max_rows=Nnodes)

        if nodes.shape[1] != 3:
            raise ValueError("Node block must have exactly 3 columns (X, Y, Z).")

        # --- Read element connectivity ---
        elements = np.loadtxt(f, max_rows=Nelem, dtype=int)

        if elements.shape[1] != 4:
            raise ValueError("Element block must have exactly 4 node IDs per element.")

    return {'nodes': nodes,
            'elements':elements}