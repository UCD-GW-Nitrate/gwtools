import pandas as pd
import numpy as np
import geopandas as gpd
from collections import defaultdict
from shapely.geometry import Polygon, Point, LineString
import plotly.graph_objects as go

def mesh_to_gdf(nodes, elements, crs):
    node_lookup = nodes.set_index("ID")[["X", "Y"]]
    geometries = []
    ies = []
    irges = []

    for _, row in elements.iterrows():
        nids = [row["ND1"], row["ND2"], row["ND3"]]
        if row["ND4"] != 0:
            nids.append(row["ND4"])

        # Convert node IDs to coordinates
        coords = [(node_lookup.loc[nid, "X"], node_lookup.loc[nid, "Y"]) for nid in nids]

        # Create polygon
        poly = Polygon(coords)

        geometries.append(poly)
        ies.append(row["IE"])
        irges.append(row["IRGE"])

    gdf = gpd.GeoDataFrame(
        {"IE": ies, "IRGE": irges},
        geometry=geometries,
        crs=crs
    )

    return gdf

def get_date_df(start_date, end_date,freq="MS"):
    c2vs_dates = pd.date_range(start=start_date, end=end_date, freq=freq)
    c2vs_dates_df = pd.DataFrame({
        'date': c2vs_dates,
        'year': c2vs_dates.year,
        'month': c2vs_dates.month,
        'days_in_month': c2vs_dates.days_in_month
    })
    return c2vs_dates_df

def stratigraphy_to_elevations(
    strat_df,
    conversion_factor=0.3048,
    id_col="ID",
    gse_col="GSE",
    remove_zero_layers=True,
    verbose=True,
):
    """
        Convert IWFM-style stratigraphy thickness columns to elevation surfaces.

        Parameters
        ----------
        strat_df : DataFrame
            Stratigraphy table with fields:
            ID, GSE, A1, L1, A2, L2, ...
        conversion_factor : float, default=0.3048
            Unit conversion factor. Default converts feet to meters.
            Use 1.0 if the input is already in desired units.
        id_col : str, default="ID"
            Node/point ID field.
        gse_col : str, default="GSE"
            Ground surface elevation field.
        remove_zero_layers : bool, default=True
            If True, remove elevation surfaces where the corresponding
            layer thickness is zero for all rows.
        verbose : bool, default=True
            If True, print reminder messages.

        Returns
        -------
        ids : ndarray
            Original ID values.
        elev : ndarray
            Elevation matrix with shape (n_points, n_surfaces).
            First column is GSE.
        elev_cols : list[str]
            Names of elevation surfaces.
        """

    ids = strat_df[id_col].to_numpy()

    # Detect A/L columns
    a_cols = sorted(
        [c for c in strat_df.columns if c.startswith("A")],
        key=lambda x: int(x[1:])
    )
    l_cols = sorted(
        [c for c in strat_df.columns if c.startswith("L")],
        key=lambda x: int(x[1:])
    )

    thickness_cols = []
    for a, l in zip(a_cols, l_cols):
        thickness_cols.extend([a, l])

    # Thickness matrix
    thickness = strat_df[thickness_cols].to_numpy(dtype=float)

    # Convert units
    gse = strat_df[gse_col].to_numpy(dtype=float) * conversion_factor
    thickness = thickness * conversion_factor

    # Build elevation surfaces: GSE, GSE-A1, GSE-A1-L1, ...
    elev = np.column_stack([
        gse,
        gse[:, None] - np.cumsum(thickness, axis=1)
    ])

    elev_cols = [gse_col] + [f"BOT_{c}" for c in thickness_cols]

    # Remove surfaces associated with layers that are zero everywhere
    if remove_zero_layers:
        keep = [True]  # always keep GSE

        for j in range(thickness.shape[1]):
            keep.append(not np.all(thickness[:, j] == 0))

        keep = np.array(keep, dtype=bool)

        elev = elev[:, keep]
        elev_cols = [c for c, k in zip(elev_cols, keep) if k]

    # Replace zero-thickness elevations by average of above and below
    corrected_count = 0

    for j in range(1, elev.shape[1] - 1):
        zero_mask = np.isclose(elev[:, j], elev[:, j - 1])

        if np.any(zero_mask):
            elev[zero_mask, j] = 0.5 * (
                    elev[zero_mask, j - 1] +
                    elev[zero_mask, j + 1]
            )
            corrected_count += int(np.sum(zero_mask))


    if verbose and corrected_count > 0:
        print(
            f"Reminder: {corrected_count} zero-thickness elevation values "
            "were replaced by the mean elevation of the surfaces above and below."
        )

    return ids, elev, elev_cols


def embed_hk_anomaly_to_nodes(
    HVK, # Horizontal vertical K properties
    HKanomaly,
    elem_df,
    id_col="ID",
    hk_col="PKH",
    elem_id_col="IE",
    elem_node_cols=("ND1", "ND2", "ND3", "ND4"),
    anomaly_elem_col=1,
    anomaly_first_layer_col=2,
):
    """
    Embed element-based HK anomaly values into node-based HKV dataframes.

    Parameters
    ----------
    HVK : list of DataFrame
        One dataframe per layer. Each dataframe must contain node IDs and
        hydraulic conductivity values, e.g. ID and PKH.
    HKanomaly : ndarray
        Array where:
        column 0 = unused ID,
        column 1 = element ID,
        columns 2 onward = anomaly values per layer.
    elem_df : DataFrame
        Element connectivity table with fields IE, ND1, ND2, ND3, ND4.
        If ND4 == 0, the element is treated as a triangle.

    Returns
    -------
    HKV_new : list of DataFrame
        Updated copy of HKV with anomaly values assigned to nodes touching
        anomalous elements.
    """

    nlay = len(HVK)

    if HKanomaly.shape[1] - anomaly_first_layer_col != nlay:
        raise ValueError(
            "Number of anomaly layer columns does not match len(HKV). "
            f"Got {HKanomaly.shape[1] - anomaly_first_layer_col} anomaly layers "
            f"and {nlay} HKV layers."
        )

    elem_node_cols = [c for c in elem_node_cols if c in elem_df.columns]

    # Element ID -> node IDs
    elem_conn = elem_df.set_index(elem_id_col)[elem_node_cols]

    # Work on copies so original HKV is not modified
    HKV_new = [df.copy() for df in HVK]

    # For fast node lookup in each layer dataframe
    node_index_maps = [
        pd.Series(df.index.to_numpy(), index=df[id_col]).to_dict()
        for df in HKV_new
    ]

    for row in HKanomaly:
        ie = row[anomaly_elem_col]

        if ie not in elem_conn.index:
            continue

        nodes = elem_conn.loc[ie].to_numpy()
        nodes = nodes[(~pd.isna(nodes)) & (nodes != 0)]

        for ilay in range(nlay):
            anomaly_value = row[anomaly_first_layer_col + ilay]

            idx_map = node_index_maps[ilay]

            for node_id in nodes:
                if node_id in idx_map:
                    HKV_new[ilay].at[idx_map[node_id], hk_col] = anomaly_value

    return HKV_new


def hvk_to_arrays(HVK, props):
    """
    Convert a layer-wise HVK list into node x layer numpy arrays.

    Parameters
    ----------
    HVK : list of DataFrame
        One dataframe per layer.
    props : str or list[str]
        Property name(s) to extract.

    Returns
    -------
    ndarray or tuple(ndarray)
        If props is a string:
            arr : (nnodes, nlayers)

        If props is a list:
            tuple of arrays in the same order as props.
    """

    single_prop = isinstance(props, str)

    if single_prop:
        props = [props]

    nlay = len(HVK)

    # sanity checks
    for p in props:
        if p not in HVK[0].columns:
            raise KeyError(f"Property '{p}' not found in HVK.")

    nnodes = len(HVK[0])

    out = {}

    for prop in props:
        arr = np.empty((nnodes, nlay), dtype=float)

        for ilay, df in enumerate(HVK):
            arr[:, ilay] = df[prop].to_numpy()

        out[prop] = arr

    if single_prop:
        return out[props[0]]

    return tuple(out[p] for p in props)


def watersheds_to_gdfs(watersheds, nodes_df, crs):
    """
    Convert IWFM small watershed connectivity to GeoDataFrames.

    Parameters
    ----------
    watersheds : list[dict]
        Output from read_iwfm_smallwatersheds().

    nodes_df : pandas.DataFrame
        DataFrame with fields:
        ID, X, Y

    crs : str, dict, pyproj.CRS, or int
        CRS passed directly to GeoDataFrame.

    Returns
    -------
    line_gdf : GeoDataFrame
    point_gdf : GeoDataFrame
    """

    required_cols = {"ID", "X", "Y"}
    missing = required_cols - set(nodes_df.columns)
    if missing:
        raise ValueError(f"nodes_df is missing required columns: {missing}")

    node_lookup = nodes_df.set_index("ID")[["X", "Y"]].to_dict("index")

    line_records = []
    point_records = []

    for ws in watersheds:
        ws_id = ws["ID"]
        iwb_list = ws["IWB"]
        qmaxwb_list = ws["QMAXWB"]

        if len(iwb_list) != len(qmaxwb_list):
            raise ValueError(
                f"Watershed ID {ws_id}: IWB and QMAXWB have different lengths"
            )

        base_attrs = {
            key: value
            for key, value in ws.items()
            if key not in {"IWB", "QMAXWB"}
        }

        coords = []

        for iwb, qmaxwb in zip(iwb_list, qmaxwb_list):
            if iwb not in node_lookup:
                raise KeyError(
                    f"Watershed ID {ws_id}: IWB node {iwb} not found in nodes_df"
                )

            x = node_lookup[iwb]["X"]
            y = node_lookup[iwb]["Y"]

            coords.append((x, y))

            point_records.append({
                **base_attrs,
                "IWB": iwb,
                "QMAXWB": qmaxwb,
                "geometry": Point(x, y),
            })

        # Skip watersheds with only one IWB
        if len(coords) > 1:
            line_records.append({
                **base_attrs,
                "IWB": iwb_list,
                "geometry": LineString(coords),
            })

    line_gdf = gpd.GeoDataFrame(line_records, geometry="geometry",crs=crs)
    point_gdf = gpd.GeoDataFrame(point_records, geometry="geometry",crs=crs)

    return line_gdf, point_gdf


def get_mesh_boundary(node, elem, crs="EPSG:26910",
                      boundary_type="outer"):
    """
    Extract ordered boundary nodes and edges from an IWFM mesh.

    Supports triangular elements (ND4=0), quadrilateral elements,
    and internal holes.

    Parameters
    ----------
    node : DataFrame
        Columns: ID, X, Y

    elem : DataFrame
        Columns: IE, ND1, ND2, ND3, ND4

    crs : str or int
        Coordinate reference system.

    boundary_type : {"outer", "all", "holes"}
        "outer" : return the outer perimeter (default).
        "all"   : return lists of all boundary loops.
        "holes" : return lists of internal boundary loops.

    Returns
    -------
    For boundary_type="outer":
        boundary_nodes, boundary_lines

    For boundary_type="all" or "holes":
        node_loops, line_loops

    boundary_nodes fields:
        ID        : global node ID
        X, Y      : coordinates
        local_id  : zero-based boundary node index
        distance  : cumulative distance along boundary
        geometry  : Point

    boundary_lines fields:
        ID        : global element ID
        nd_id0    : global starting node ID
        nd_id1    : global ending node ID
        nd_loc0   : local starting node index
        nd_loc1   : local ending node index
        length    : edge length
        distance  : cumulative distance at edge start
        geometry  : LineString

    Notes
    -----
    Assumes a connected mesh domain with one outer boundary
    and zero or more holes.

    All boundary loops must be simple and closed.
    """

    if boundary_type not in ("outer", "all", "holes"):
        raise ValueError(
            "boundary_type must be 'outer', 'all', or 'holes'"
        )

    # -------------------------------------------------------
    # 1. Extract triangle and quadrilateral edges
    # -------------------------------------------------------
    elem_nodes = elem[
        ["ND1", "ND2", "ND3", "ND4"]
    ].to_numpy(dtype=int)

    elem_ids = elem["IE"].to_numpy(dtype=int)

    is_triangle = elem_nodes[:, 3] == 0
    is_quad = ~is_triangle

    tri_pairs = np.array([
        [0, 1],
        [1, 2],
        [2, 0],
    ])

    quad_pairs = np.array([
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 0],
    ])

    tri_edges = elem_nodes[is_triangle][:, tri_pairs].reshape(-1, 2)
    quad_edges = elem_nodes[is_quad][:, quad_pairs].reshape(-1, 2)

    edges = np.vstack([tri_edges, quad_edges])

    edge_elem_ids = np.concatenate([
        np.repeat(elem_ids[is_triangle], 3),
        np.repeat(elem_ids[is_quad], 4),
    ])

    if np.any(edges <= 0):
        raise ValueError("Invalid node ID found in mesh edges.")

    # -------------------------------------------------------
    # 2. Identify boundary edges
    # -------------------------------------------------------
    sorted_edges = np.sort(edges, axis=1)

    _, inverse, counts = np.unique(
        sorted_edges,
        axis=0,
        return_inverse=True,
        return_counts=True,
    )

    if np.any(counts > 2):
        raise ValueError(
            "Non-manifold mesh: some edges belong to more "
            "than two elements."
        )

    boundary_mask = counts[inverse] == 1

    boundary_edges = edges[boundary_mask]
    boundary_elem_ids = edge_elem_ids[boundary_mask]

    if len(boundary_edges) == 0:
        raise ValueError("No boundary edges found.")

    # -------------------------------------------------------
    # 3. Build adjacency
    # -------------------------------------------------------
    adjacency = defaultdict(list)

    for i, (n0, n1) in enumerate(boundary_edges):
        adjacency[int(n0)].append((int(n1), i))
        adjacency[int(n1)].append((int(n0), i))

    invalid = [
        n for n, neighbors in adjacency.items()
        if len(neighbors) != 2
    ]

    if invalid:
        raise ValueError(
            f"Boundary contains {len(invalid)} nodes "
            f"with degree != 2. Examples: {invalid[:10]}"
        )

    # -------------------------------------------------------
    # 4. Traverse all closed boundary loops
    # -------------------------------------------------------
    visited_edges = set()
    loops = []

    for start_edge in range(len(boundary_edges)):

        if start_edge in visited_edges:
            continue

        start_node = int(boundary_edges[start_edge, 0])
        current = start_node
        previous = None

        ordered_nodes = []
        ordered_edges = []

        while True:

            ordered_nodes.append(current)

            candidates = [
                (n, idx)
                for n, idx in adjacency[current]
                if idx not in visited_edges
            ]

            if not candidates:
                raise ValueError(
                    "Boundary traversal failed: open or "
                    "inconsistent boundary."
                )

            if previous is None:
                next_node, edge_idx = next(
                    ((n, idx) for n, idx in candidates
                     if idx == start_edge),
                    candidates[0]
                )
            else:
                next_node, edge_idx = next(
                    ((n, idx) for n, idx in candidates
                     if n != previous),
                    candidates[0]
                )

            visited_edges.add(edge_idx)
            ordered_edges.append(edge_idx)

            if next_node == start_node:
                break

            previous = current
            current = next_node

        loops.append((
            np.asarray(ordered_nodes, dtype=int),
            np.asarray(ordered_edges, dtype=int)
        ))

    # -------------------------------------------------------
    # 5. Determine outer loop by enclosed area
    # -------------------------------------------------------
    node_lookup = node.set_index("ID")

    loop_areas = []

    for ordered_nodes, _ in loops:

        xy = node_lookup.loc[
            ordered_nodes, ["X", "Y"]
        ].to_numpy()

        poly = Polygon(xy)

        if not poly.is_valid:
            raise ValueError(
                "Invalid boundary polygon encountered."
            )

        loop_areas.append(poly.area)

    outer_idx = int(np.argmax(loop_areas))

    # Outer loop first, remaining loops by decreasing area
    loop_order = [outer_idx] + sorted(
        [i for i in range(len(loops)) if i != outer_idx],
        key=lambda i: loop_areas[i],
        reverse=True
    )

    # -------------------------------------------------------
    # 6. Build ordered GeoDataFrames
    # -------------------------------------------------------
    node_loops = []
    line_loops = []

    for loop_idx in loop_order:

        ordered_nodes, ordered_edges = loops[loop_idx]

        n = len(ordered_nodes)

        # Ordered boundary nodes
        bnodes = node_lookup.loc[
            ordered_nodes, ["X", "Y"]
        ].copy()

        bnodes.insert(0, "ID", ordered_nodes)

        bnodes["local_id"] = np.arange(n, dtype=int)

        bnodes = gpd.GeoDataFrame(
            bnodes.reset_index(drop=True),
            geometry=gpd.points_from_xy(
                bnodes["X"], bnodes["Y"]
            ),
            crs=crs,
        )

        # Ordered boundary edges
        xy = bnodes[["X", "Y"]].to_numpy()
        xy_next = np.roll(xy, -1, axis=0)

        lengths = np.linalg.norm(
            xy_next - xy, axis=1
        )

        distances = np.concatenate([
            [0.0],
            np.cumsum(lengths[:-1])
        ])

        bnodes["distance"] = distances

        blines = gpd.GeoDataFrame(
            {
                "ID": boundary_elem_ids[ordered_edges],
                "nd_id0": ordered_nodes,
                "nd_id1": np.roll(ordered_nodes, -1),
                "nd_loc0": np.arange(n, dtype=int),
                "nd_loc1": np.roll(
                    np.arange(n, dtype=int), -1
                ),
                "length": lengths,
                "distance": distances,
            },
            geometry=[
                LineString([p0, p1])
                for p0, p1 in zip(xy, xy_next)
            ],
            crs=crs,
        )

        node_loops.append(bnodes)
        line_loops.append(blines)

    # -------------------------------------------------------
    # 7. Return requested boundary type
    # -------------------------------------------------------
    if boundary_type == "outer":
        return node_loops[0], line_loops[0]

    if boundary_type == "holes":
        return node_loops[1:], line_loops[1:]

    return node_loops, line_loops


def get_boundary_parts(bnd_nodes, bnd_lines, field, threshold,
                       use_abs=None):
    """
    Identify continuous boundary parts satisfying a field threshold.

    Parameters
    ----------
    bnd_nodes : GeoDataFrame
        Ordered boundary nodes with the selected statistic.
        Must contain 'ID' and 'local_id'.

    bnd_lines : GeoDataFrame
        Ordered boundary edges with 'nd_loc0', 'nd_loc1', 'length'.

    field : str
        Field to evaluate, e.g. 'slope_monthly_L1' or 'std_yearly_L2'.

    threshold : float
        Maximum allowed value (strictly less than).

    use_abs : bool, optional
        If None, automatically True for slope fields and False
        for other fields.

    Returns
    -------
    parts : list of dict
        Each dictionary contains:
            'nodes'       : ordered node GeoDataFrame
            'lines'       : ordered line GeoDataFrame
            'length'      : total boundary length (m)
            'length_km'   : total boundary length (km)
            'n_nodes'     : number of nodes
            'n_segments'  : number of edges

    Notes
    -----
    An edge is valid only if BOTH endpoints satisfy the threshold.
    Isolated valid nodes are excluded.
    The boundary is treated as a closed loop.
    """

    if use_abs is None:
        use_abs = field.startswith("slope")

    if threshold < 0:
        raise ValueError("threshold must be non-negative")

    n = len(bnd_nodes)

    if len(bnd_lines) != n:
        raise ValueError(
            "Expected one boundary edge per boundary node."
        )

    # Verify ordered boundary topology
    idx = np.arange(n)

    if not np.array_equal(
        bnd_lines["nd_loc0"].to_numpy(), idx
    ) or not np.array_equal(
        bnd_lines["nd_loc1"].to_numpy(), (idx + 1) % n
    ):
        raise ValueError(
            "Boundary nodes and lines are not consistently ordered."
        )

    # -------------------------------------------------------
    # 1. Identify valid nodes
    # -------------------------------------------------------
    values = bnd_nodes[field].to_numpy(dtype=float)

    if use_abs:
        valid_nodes = np.isfinite(values) & (
            np.abs(values) < threshold
        )
    else:
        valid_nodes = np.isfinite(values) & (
            values < threshold
        )

    # -------------------------------------------------------
    # 2. Identify valid edges
    # -------------------------------------------------------
    # Edge i connects node i to node (i+1) % n
    valid_edges = valid_nodes & np.roll(valid_nodes, -1)

    if not valid_edges.any():
        return []

    # -------------------------------------------------------
    # 3. Group consecutive valid edges (circular boundary)
    # -------------------------------------------------------
    if valid_edges.all():
        edge_groups = [np.arange(n)]
    else:
        # A valid edge preceded by an invalid edge starts a part
        starts = np.where(
            valid_edges & ~np.roll(valid_edges, 1)
        )[0]

        edge_groups = []

        for start in starts:
            group = []
            i = int(start)

            while valid_edges[i]:
                group.append(i)
                i = (i + 1) % n

            edge_groups.append(np.asarray(group, dtype=int))

    # -------------------------------------------------------
    # 4. Build output dictionaries
    # -------------------------------------------------------
    parts = []

    for group in edge_groups:

        # Consecutive edges: collect starting nodes and final node
        if len(group) == n:
            node_idx = np.arange(n)
        else:
            node_idx = np.r_[
                group,
                (group[-1] + 1) % n
            ]

        part_nodes = bnd_nodes.iloc[node_idx].copy()
        part_lines = bnd_lines.iloc[group].copy()

        # Reset row indices, preserve original local_id
        part_nodes = part_nodes.reset_index(drop=True)
        part_lines = part_lines.reset_index(drop=True)

        # Distance measured from start of this part
        edge_lengths = part_lines["length"].to_numpy()

        if len(group) == n:
            node_dist = np.r_[
                0.0, np.cumsum(edge_lengths[:-1])
            ]
        else:
            node_dist = np.r_[
                0.0, np.cumsum(edge_lengths)
            ]

        part_nodes["part_distance"] = node_dist
        part_lines["part_distance"] = np.r_[
            0.0, np.cumsum(edge_lengths[:-1])
        ]

        total_length = float(edge_lengths.sum())

        parts.append({
            "nodes": part_nodes,
            "lines": part_lines,
            "length": total_length,
            "length_km": total_length / 1000.0,
            "n_nodes": len(part_nodes),
            "n_segments": len(part_lines),
        })

    # Sort longest boundary parts first
    parts.sort(key=lambda p: p["length"], reverse=True)

    return parts


def get_stable_boundary_parts(
    bnd_nodes,
    bnd_lines,
    slope_field,
    std_field,
    slope_threshold,
    std_threshold,
):
    """
    Identify continuous boundary parts satisfying BOTH:
        abs(slope) < slope_threshold
        std        < std_threshold

    A boundary edge is valid only when both endpoint nodes
    satisfy both conditions.

    Assumes a single closed, ordered boundary loop.

    Returns
    -------
    parts : list of dictionaries
        nodes       : GeoDataFrame of ordered nodes
        lines       : GeoDataFrame of ordered boundary edges
        length      : total length (m)
        length_km   : total length (km)
        n_nodes     : number of nodes
        n_segments  : number of edges
    """

    if slope_threshold < 0 or std_threshold < 0:
        raise ValueError("Thresholds must be non-negative.")

    n = len(bnd_nodes)

    if len(bnd_lines) != n:
        raise ValueError("Expected one edge per boundary node.")

    # Verify topology
    idx = np.arange(n)

    if not np.array_equal(
        bnd_lines["nd_loc0"].to_numpy(), idx
    ) or not np.array_equal(
        bnd_lines["nd_loc1"].to_numpy(), (idx + 1) % n
    ):
        raise ValueError(
            "Boundary nodes and lines are not consistently ordered."
        )

    # -------------------------------------------------------
    # 1. Valid nodes: both criteria must be satisfied
    # -------------------------------------------------------
    slope = bnd_nodes[slope_field].to_numpy(dtype=float)
    std = bnd_nodes[std_field].to_numpy(dtype=float)

    valid_nodes = (
        np.isfinite(slope)
        & np.isfinite(std)
        & (np.abs(slope) < slope_threshold)
        & (std >= 0)
        & (std < std_threshold)
    )

    # -------------------------------------------------------
    # 2. Valid edges: both endpoints must be valid
    # -------------------------------------------------------
    valid_edges = valid_nodes & np.roll(valid_nodes, -1)

    if not valid_edges.any():
        return []

    # -------------------------------------------------------
    # 3. Identify continuous edge groups (circular boundary)
    # -------------------------------------------------------
    if valid_edges.all():
        edge_groups = [np.arange(n)]
    else:
        starts = np.where(
            valid_edges & ~np.roll(valid_edges, 1)
        )[0]

        edge_groups = []

        for start in starts:
            group = []
            i = int(start)

            while valid_edges[i]:
                group.append(i)
                i = (i + 1) % n

            edge_groups.append(np.asarray(group, dtype=int))

    # -------------------------------------------------------
    # 4. Construct boundary parts
    # -------------------------------------------------------
    parts = []

    for group in edge_groups:

        if len(group) == n:
            node_idx = np.arange(n)
        else:
            node_idx = np.r_[
                group,
                (group[-1] + 1) % n
            ]

        part_nodes = (
            bnd_nodes.iloc[node_idx]
            .copy()
            .reset_index(drop=True)
        )

        part_lines = (
            bnd_lines.iloc[group]
            .copy()
            .reset_index(drop=True)
        )

        edge_lengths = part_lines["length"].to_numpy(dtype=float)

        if len(group) == n:
            node_dist = np.r_[
                0.0,
                np.cumsum(edge_lengths[:-1])
            ]
        else:
            node_dist = np.r_[
                0.0,
                np.cumsum(edge_lengths)
            ]

        part_nodes["part_distance"] = node_dist

        part_lines["part_distance"] = np.r_[
            0.0,
            np.cumsum(edge_lengths[:-1])
        ]

        total_length = float(edge_lengths.sum())

        parts.append({
            "nodes": part_nodes,
            "lines": part_lines,
            "length": total_length,
            "length_km": total_length / 1000.0,
            "n_nodes": len(part_nodes),
            "n_segments": len(part_lines),
        })

    parts.sort(key=lambda p: p["length"], reverse=True)

    return parts


def plot_stable_boundary_parts(
    bnd_nodes,
    bnd_lines,
    parts,
    layer=1,
    color_field=None,
    marker_size=6,
    fig=None,
):
    """
    Optimized Plotly visualization of stable boundary parts.

    Parameters
    ----------
    bnd_nodes : GeoDataFrame
        Original boundary nodes with statistics.

    bnd_lines : GeoDataFrame
        Original boundary edges.

    parts : list of dict
        Output from get_stable_boundary_parts().

    layer : int
        C2VSim layer (1-4), used for hover statistics.

    color_field : str, optional
        Field used to color stable boundary nodes.
        Default: std_yearly_L{layer}.

    marker_size : float
        Marker size.

    fig : go.Figure, optional
        Existing Plotly figure.

    Returns
    -------
    fig : go.Figure
    """

    std_field = f"std_yearly_L{layer}"
    slope_field = f"slope_yearly_L{layer}"

    if color_field is None:
        color_field = std_field

    if fig is None:
        fig = go.Figure()

    # -------------------------------------------------------
    # 1. Entire CV boundary (single trace)
    # -------------------------------------------------------
    x_all, y_all = [], []

    for line in bnd_lines.geometry:
        x, y = line.xy

        x_all.extend(x)
        x_all.append(None)

        y_all.extend(y)
        y_all.append(None)

    fig.add_trace(go.Scatter(
        x=x_all,
        y=y_all,
        mode="lines",
        line=dict(color="lightgray", width=1),
        hoverinfo="skip",
        name="CV boundary",
        showlegend=True,
    ))

    # -------------------------------------------------------
    # 2. Stable boundary lines (single trace)
    # -------------------------------------------------------
    x_stable, y_stable = [], []

    for part in parts:
        for line in part["lines"].geometry:
            x, y = line.xy

            x_stable.extend(x)
            x_stable.append(None)

            y_stable.extend(y)
            y_stable.append(None)

    if x_stable:
        fig.add_trace(go.Scatter(
            x=x_stable,
            y=y_stable,
            mode="lines",
            line=dict(color="black", width=3),
            hoverinfo="skip",
            name="Stable boundary",
            showlegend=True,
        ))

    # -------------------------------------------------------
    # 3. Combine all stable nodes into one dataset
    # -------------------------------------------------------
    if parts:

        node_frames = []

        for i, part in enumerate(parts):

            tmp = part["nodes"].copy()

            tmp["part_id"] = i + 1
            tmp["part_length_km"] = part["length_km"]

            node_frames.append(tmp)

        import pandas as pd

        stable_nodes = pd.concat(
            node_frames,
            ignore_index=True
        )

        # ---------------------------------------------------
        # 4. Color scale
        # ---------------------------------------------------
        values = stable_nodes[color_field].to_numpy(dtype=float)

        finite = values[np.isfinite(values)]

        if len(finite):
            cmin = float(finite.min())
            cmax = float(finite.max())

            if color_field.startswith("slope"):
                vmax = max(abs(cmin), abs(cmax))
                cmin, cmax = -vmax, vmax
                colorscale = "RdBu_r"
            else:
                colorscale = "Viridis"

            if cmin == cmax:
                cmin -= 0.5
                cmax += 0.5
        else:
            cmin, cmax = 0, 1
            colorscale = "Viridis"

        # ---------------------------------------------------
        # 5. Hover information
        # ---------------------------------------------------
        customdata = np.column_stack([
            stable_nodes["part_id"].to_numpy(),
            stable_nodes["ID"].to_numpy(),
            stable_nodes["distance"].to_numpy() / 1000,
            stable_nodes["part_distance"].to_numpy() / 1000,
            stable_nodes["part_length_km"].to_numpy(),
            stable_nodes[std_field].to_numpy(),
            stable_nodes[slope_field].to_numpy(),
        ])

        hovertemplate = (
            "<b>Part %{customdata[0]:.0f}</b><br>"
            "Node ID: %{customdata[1]:.0f}<br>"
            "Boundary distance: %{customdata[2]:.2f} km<br>"
            "Part distance: %{customdata[3]:.2f} km<br>"
            "Part length: %{customdata[4]:.2f} km<br><br>"
            f"<b>Layer {layer}</b><br>"
            "STD yearly: %{customdata[5]:.3f} m<br>"
            "Slope yearly: %{customdata[6]:+.3f} m/year<br>"
            f"{color_field}: %{{marker.color:.3f}}"
            "<extra></extra>"
        )

        # ---------------------------------------------------
        # 6. Stable nodes (single trace)
        # ---------------------------------------------------
        fig.add_trace(go.Scatter(
            x=stable_nodes.geometry.x.to_numpy(),
            y=stable_nodes.geometry.y.to_numpy(),
            mode="markers",
            customdata=customdata,
            marker=dict(
                size=marker_size,
                color=values,
                colorscale=colorscale,
                cmin=cmin,
                cmax=cmax,
                showscale=True,
                colorbar=dict(title=color_field),
                line=dict(width=0),
            ),
            hovertemplate=hovertemplate,
            name="Stable nodes",
            showlegend=False,
        ))

    # -------------------------------------------------------
    # 7. Summary
    # -------------------------------------------------------
    total_length = sum(p["length"] for p in parts)
    cv_length = bnd_lines["length"].sum()

    percentage = (
        100 * total_length / cv_length
        if cv_length > 0 else 0
    )

    fig.update_layout(
        title=(
            f"Stable boundary parts — Layer {layer}"
            f"<br><sup>{len(parts)} parts | "
            f"{total_length / 1000:.1f} km | "
            f"{percentage:.1f}% of CV boundary</sup>"
        ),
        xaxis_title="Easting (m)",
        yaxis_title="Northing (m)",
        #template="plotly_white",
        width=1100,
        height=850,
        hovermode="closest",
        legend=dict(
            x=0.01,
            y=0.99,
            xanchor="left",
            yanchor="top",
        ),
    )

    fig.update_yaxes(
        scaleanchor="x",
        scaleratio=1
    )

    return fig


def assign_boundary_head_elevations(
    edited_parts,
    bnd_nodes_stats,
    strat_elev,
):
    """
    Assign global row IDs to edited boundary nodes and build
    separate arrays for NPSAT boundary heads and elevations.

    Parameters
    ----------
    edited_parts : list of dict
        Edited projected boundary parts.
        Each part contains a 'nodes' GeoDataFrame with original
        C2VSim node IDs in column 'ID'.

    bnd_nodes_stats : GeoDataFrame
        Original boundary statistics containing:
            ID, mean_L1, mean_L2, mean_L3, mean_L4

    strat_elev : tuple
        Output of iwfm_u.stratigraphy_to_elevations():
            ids, elev, elev_cols

    Returns
    -------
    parts_out : list of dict
        Copies of edited_parts with 'row_id' added to nodes.

    bc_heads : ndarray, shape (N, 4)
        Mean heads for the four C2VSim layers (m).

    bc_elev : ndarray, shape (N, 3)
        Interlayer elevations:
            BOT_L1, BOT_L2, BOT_L3 (m).

    Notes
    -----
    row_id is a global zero-based index across all parts.

    Nodes appearing in multiple parts receive separate row IDs.
    No spatial interpolation is performed.
    """

    # -------------------------------------------------------
    # 1. Extract stratigraphy arrays
    # -------------------------------------------------------
    strat_ids, elevations, elev_cols = strat_elev

    strat_ids = np.asarray(strat_ids, dtype=int)
    elevations = np.asarray(elevations, dtype=float)

    head_cols = [
        "mean_L1",
        "mean_L2",
        "mean_L3",
        "mean_L4",
    ]

    required_elev_cols = [
        "BOT_L1",
        "BOT_L2",
        "BOT_L3",
    ]

    elev_indices = [
        elev_cols.index(col)
        for col in required_elev_cols
    ]

    # -------------------------------------------------------
    # 2. Create lookup tables using original node IDs
    # -------------------------------------------------------
    stats_ids = bnd_nodes_stats["ID"].to_numpy(dtype=int)

    if len(np.unique(stats_ids)) != len(stats_ids):
        raise ValueError("Duplicate node IDs in bnd_nodes_stats.")

    if len(np.unique(strat_ids)) != len(strat_ids):
        raise ValueError("Duplicate node IDs in strat_elev.")

    if elevations.shape[0] != len(strat_ids):
        raise ValueError("Stratigraphy ID/elevation shape mismatch.")

    stats_lookup = {
        int(node_id): i
        for i, node_id in enumerate(stats_ids)
    }

    strat_lookup = {
        int(node_id): i
        for i, node_id in enumerate(strat_ids)
    }

    head_values = bnd_nodes_stats[
        head_cols
    ].to_numpy(dtype=float)

    elev_values = elevations[:, elev_indices]

    # -------------------------------------------------------
    # 3. Assign global row IDs and build arrays
    # -------------------------------------------------------
    parts_out = []

    heads_list = []
    elev_list = []

    row_offset = 0

    for part in edited_parts:

        new_part = part.copy()
        nodes = part["nodes"].copy().reset_index(drop=True)

        node_ids = nodes["ID"].to_numpy(dtype=int)

        missing_stats = [
            nid for nid in node_ids
            if nid not in stats_lookup
        ]

        missing_strat = [
            nid for nid in node_ids
            if nid not in strat_lookup
        ]

        if missing_stats or missing_strat:
            raise ValueError(
                f"Missing node IDs. "
                f"Stats: {missing_stats[:10]}, "
                f"Stratigraphy: {missing_strat[:10]}"
            )

        stats_idx = np.array([
            stats_lookup[nid] for nid in node_ids
        ], dtype=int)

        strat_idx = np.array([
            strat_lookup[nid] for nid in node_ids
        ], dtype=int)

        n = len(nodes)

        nodes["row_id"] = np.arange(
            row_offset,
            row_offset + n,
            dtype=int
        )

        heads_list.append(head_values[stats_idx])
        elev_list.append(elev_values[strat_idx])

        new_part["nodes"] = nodes
        parts_out.append(new_part)

        row_offset += n

    # -------------------------------------------------------
    # 4. Assemble final arrays
    # -------------------------------------------------------
    bc_heads = (
        np.vstack(heads_list)
        if heads_list
        else np.empty((0, 4), dtype=float)
    )

    bc_elev = (
        np.vstack(elev_list)
        if elev_list
        else np.empty((0, 3), dtype=float)
    )

    if not np.isfinite(bc_heads).all():
        raise ValueError("Non-finite values found in bc_heads.")

    if not np.isfinite(bc_elev).all():
        raise ValueError("Non-finite values found in bc_elev.")

    return parts_out, bc_heads, bc_elev