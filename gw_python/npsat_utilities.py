import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon, LineString

def read_urf_calc(prefix_filename, nproc, n_strml, por, show_progress=False):
    col_names = ['SrcID', 'SrcInd', 'Eid', 'Sid', 'ER',
                 'p_cdsX', 'p_cdsY', 'p_cdsZ', 'v_cds',
                 'p_lndX', 'p_lndY', 'Len']

    urf_data = pd.DataFrame(np.nan, index=range(n_strml), columns=col_names)
    urf_msa = np.full((n_strml, 3 * len(por)), np.nan)

    ida = 0
    for ii in range(nproc):
        fname = prefix_filename + f"{ii}.dat"
        #print("update2")
        if show_progress:
            print(f"Reading {fname}")
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
    streams = []

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
            coords = np.array(
                [list(map(float, f.readline().strip().split()))
                 for _ in range(nv)],
                dtype=float,
            )

            # Store original coordinates
            streams.append(coords)

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

    return gdf, streams

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

def polygon_to_avg_line(gdf):
    """
    From a polygon GeoDataFrame, create a new line GeoDataFrame where:
    - First point = average of vertices 0 and 3
    - Second point = average of vertices 1 and 2
    """

    new_geoms = []

    for idx, row in gdf.iterrows():
        poly = row.geometry
        exterior = list(poly.exterior.coords)

        if len(exterior) < 4:
            raise ValueError(f"Polygon at index {idx} has fewer than 4 vertices.")

        # Extract first 4 vertices
        p0 = exterior[0]
        p1 = exterior[1]
        p2 = exterior[2]
        p3 = exterior[3]

        # Compute averages
        avg03 = ((p0[0] + p3[0]) / 2, (p0[1] + p3[1]) / 2)
        avg12 = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)

        # Create line
        line = LineString([avg03, avg12])
        new_geoms.append(line)

    # Build new GeoDataFrame
    return gpd.GeoDataFrame(geometry=new_geoms, crs=gdf.crs)


def read_scattered(filename):
    """
    Read an ASCII "SCATTERED" data file.

    Expected format
    ----------------
    Line 1 : "SCATTERED"                      (literal keyword, mandatory)
    Line 2 : "3D" or "2D"                     (dimensionality flag)
    Line 3 : interpolation method(s)
               - 3D -> "<xy_interp> <z_interp>"   e.g. "LINEAR NEAREST"
               - 2D -> "<xy_interp>"              e.g. "LINEAR"
    Line 4 : Nnodes  Ndata  Ntriangles
    Next Nnodes lines (node/data block):
               - X  Y  V1 ELEV1  V2 ELEV2 ... Vk ELEVk [V(k+1)]
                 Ndata (from line 4) is the TOTAL number of columns after
                 X,Y -- not a pair count. Those Ndata columns are read as
                 (V, ELEV) pairs:
                   - Ndata even -> Ndata/2 complete pairs, no trailing V
                   - Ndata odd  -> (Ndata-1)/2 pairs, plus one trailing
                                   unpaired V(k+1) with no elevation
                 2D files are the special case Ndata=1: 0 pairs + 1 V,
                 i.e. plain "X Y V".
    Next Ntriangles lines : i  j  k            (1-based or 0-based node
                                                 indices forming a triangle)

    Returns
    -------
    dict with keys:
        'XY'        : (Nnodes, 2)   ndarray of X,Y coordinates
        'V'         : (Nnodes, n_value_cols) ndarray of data values, where
                      n_value_cols = Ndata//2 (+1 if Ndata is odd)
        'ELEV'      : (Nnodes, Ndata//2) ndarray of elevations
                      (empty ndarray for 2D files)
        'XYinterp'  : str, the horizontal interpolation method
        'Zinterp'   : str or None, the vertical interpolation method
                      (None for 2D files)
        'triangles' : (Ntriangles, 3) ndarray of int, triangle connectivity
    """
    with open(filename, 'r') as f:
        # keep only non-blank lines, strip trailing/leading whitespace
        raw_lines = [ln.strip() for ln in f if ln.strip() != '']

    if len(raw_lines) < 4:
        raise ValueError("File is too short to be a valid SCATTERED file.")

    # --- header -----------------------------------------------------
    keyword = raw_lines[0].split()[0].upper()
    if keyword != 'SCATTERED':
        raise ValueError(
            f"Expected first line to be 'SCATTERED', got '{raw_lines[0]}'"
        )

    dim = raw_lines[1].split()[0].upper()
    if dim not in ('2D', '3D'):
        raise ValueError(f"Expected '2D' or '3D' on line 2, got '{raw_lines[1]}'")

    interp_tokens = raw_lines[2].split()
    if dim == '3D':
        if len(interp_tokens) < 2:
            raise ValueError(
                "3D files require two interpolation methods "
                "(xy_interp z_interp) on line 3."
            )
        xy_interp, z_interp = interp_tokens[0], interp_tokens[1]
    else:
        if len(interp_tokens) < 1:
            raise ValueError("2D files require one interpolation method on line 3.")
        xy_interp, z_interp = interp_tokens[0], None

    counts = raw_lines[3].split()
    if len(counts) < 3:
        raise ValueError(
            "Line 4 must contain Nnodes, Ndata and Ntriangles."
        )
    Nnodes, Ndata, Ntriangles = (int(v) for v in counts[:3])

    # --- node / data block -------------------------------------------
    data_start = 4
    data_end = data_start + Nnodes
    node_lines = raw_lines[data_start:data_end]

    if len(node_lines) != Nnodes:
        raise ValueError(
            f"Expected {Nnodes} node data lines, found {len(node_lines)}."
        )

    XY = np.empty((Nnodes, 2), dtype=float)

    # Ndata is the TOTAL number of columns after X,Y (not a pair count).
    # Those columns follow the pattern V1 ELEV1 V2 ELEV2 ... Vk ELEVk [V(k+1)]:
    #   - Ndata even -> Ndata/2 complete (V, ELEV) pairs, no trailing V
    #   - Ndata odd  -> (Ndata-1)/2 complete pairs, plus one trailing
    #                   unpaired V(k+1) with no matching elevation
    # (2D files naturally fall out of this with Ndata=1: 0 pairs + 1 V.)
    n_pairs = Ndata // 2
    has_extra = (Ndata % 2) == 1
    n_value_cols = n_pairs + (1 if has_extra else 0)

    V = np.empty((Nnodes, n_value_cols), dtype=float)
    ELEV = np.empty((Nnodes, n_pairs), dtype=float) if dim == '3D' else np.empty((0, 0))

    for i, line in enumerate(node_lines):
        vals = [float(x) for x in line.split()]

        if len(vals) < 2:
            raise ValueError(f"Node line {i + 1} does not contain X, Y values.")

        XY[i, 0], XY[i, 1] = vals[0], vals[1]
        rest = vals[2:]

        if len(rest) != Ndata:
            raise ValueError(
                f"Node line {i + 1} expected {Ndata} entries after X,Y "
                f"(per header Ndata), found {len(rest)}."
            )

        pairs = rest[:2 * n_pairs]
        V[i, :n_pairs] = pairs[0::2]
        if dim == '3D':
            ELEV[i, :] = pairs[1::2]
        if has_extra:
            V[i, n_pairs] = rest[2 * n_pairs]  # trailing V(k+1), no elevation

    # --- triangle block -------------------------------------------------
    tri_start = data_end
    tri_end = tri_start + Ntriangles
    tri_lines = raw_lines[tri_start:tri_end]

    if len(tri_lines) != Ntriangles:
        raise ValueError(
            f"Expected {Ntriangles} triangle lines, found {len(tri_lines)}."
        )

    triangles = np.empty((Ntriangles, 3), dtype=int)
    for i, line in enumerate(tri_lines):
        vals = [int(float(x)) for x in line.split()[:3]]
        if len(vals) != 3:
            raise ValueError(f"Triangle line {i + 1} does not contain 3 indices.")
        triangles[i, :] = vals

    return {
        'XY': XY,
        'V': V,
        'ELEV': ELEV,
        'XYinterp': xy_interp,
        'Zinterp': z_interp,
        'triangles': triangles,
    }


import pandas as pd


def read_wells(filename):
    """
    Read an ASCII wells file.

    Expected format
    ----------------
    Line 1 : Nwells               (number of wells to follow; not strictly
                                    needed since the rest of the file is
                                    read directly, but validated if present)
    Remaining lines : X  Y  Top  Bottom  Q      (one well per line)

    Parameters
    ----------
    filename : str
        Path to the wells file.

    Returns
    -------
    pandas.DataFrame
        Columns: X, Y, Top, Bottom, Q
    """
    with open(filename, 'r') as f:
        lines = [ln.strip() for ln in f if ln.strip() != '']

    if len(lines) < 1:
        raise ValueError("File is empty.")

    Nwells = int(float(lines[0]))
    well_lines = lines[1:]

    if len(well_lines) != Nwells:
        raise ValueError(
            f"Header declares {Nwells} wells but found {len(well_lines)} "
            f"data lines."
        )

    rows = []
    for i, line in enumerate(well_lines):
        vals = [float(x) for x in line.split()]
        if len(vals) != 5:
            raise ValueError(
                f"Well line {i + 1} expected 5 columns (X, Y, Top, Bottom, "
                f"Q), found {len(vals)}."
            )
        rows.append(vals)

    df = pd.DataFrame(rows, columns=['X', 'Y', 'Top', 'Bottom', 'Q'])
    return df

def read_bc(filename):
    """
    Read boundary conditions from a NPSAT boundary file.

    Returns
    -------
    list of dict
        Each dictionary contains:
        {
            "type": str,
            "nv": int,
            "val": float or str,
            "vertices": np.ndarray of shape (nv, 2)
        }
    """

    boundaries = []

    with open(filename, "r") as f:
        # Number of boundary conditions
        N = int(f.readline().strip())

        for _ in range(N):

            # Boundary header: TYPE nv val
            line = f.readline().strip().split()

            bc_type = line[0]
            nv = int(line[1])

            # val may be either numeric or a string
            try:
                val = float(line[2])
            except ValueError:
                val = line[2]

            # Vertices
            vertices = np.empty((nv, 2), dtype=float)

            for i in range(nv):
                vertices[i, :] = [
                    float(v) for v in f.readline().split()
                ]

            boundaries.append({
                "type": bc_type,
                "nv": nv,
                "val": val,
                "vertices": vertices
            })

    return boundaries


def read_multi_rect(filename):
    """
    Read a MULTIRECT boundary/area file with the format:

        MULTIRECT
        N
        interp_type interp_file
        xmin ymin xmax ymax
        interp_type interp_file
        xmin ymin xmax ymax
        ...  (repeated N times)

    Returns
    -------
    list of dict, each with keys:
        'type' : str            - interpolation type (e.g. "GRIDDED")
        'file' : str            - interpolation data filename
        'bbox' : np.ndarray     - 2x2 array [[xmin, ymin], [xmax, ymax]]
    """
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    idx = 0
    keyword = lines[idx]; idx += 1
    if keyword != 'MULTIRECT':
        raise ValueError(f"Expected 'MULTIRECT' keyword, got '{keyword}'")

    n = int(lines[idx]); idx += 1

    area_list = []
    for _ in range(n):
        interp_type, interp_file = lines[idx].split(); idx += 1
        xmin, ymin, xmax, ymax = map(float, lines[idx].split()); idx += 1

        area_list.append({
            'type': interp_type,
            'file': interp_file,
            'bbox': np.array([[xmin, ymin], [xmax, ymax]])
        })

    return area_list
