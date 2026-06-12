import pandas as pd
import numpy as np
from shapely.geometry import Polygon
from shapely import affinity
import geopandas as gpd
from collections import defaultdict

def minimum_integer_angle_rectangle(geom, angle_step=1):
    """
    Find an oriented orthogonal rectangle enclosing a geometry, using integer-like
    angle search, and return rectangle coordinates.

    Parameters
    ----------
    geom : shapely geometry or GeoSeries/GeoDataFrame
        Geometry to enclose.
    angle_step : float, default 1
        Angle increment in degrees. Use 1 for integer angles.
    return_angle : bool, default False
        If True, also return the selected angle and rectangle polygon.

    Returns
    -------
    xy : ndarray, shape (5, 2)
        Rectangle exterior coordinates. Last point repeats first point.
    angle : float, optional
        Selected rotation angle in degrees.
    rect_poly : shapely Polygon, optional
        Rectangle polygon.
    """

    # Accept GeoSeries / GeoDataFrame / single shapely geometry
    if hasattr(geom, "unary_union"):
        geom = geom.unary_union

    if geom.is_empty:
        raise ValueError("Input geometry is empty.")

    best_area = np.inf
    best_rect = None
    best_angle = None

    # Search only 0-90 degrees because rectangles repeat by 90 degrees
    angles = np.arange(0, 90, angle_step)

    for angle in angles:
        # Rotate geometry so candidate rectangle is axis-aligned
        g_rot = affinity.rotate(geom, -angle, origin="centroid", use_radians=False)

        minx, miny, maxx, maxy = g_rot.bounds
        area = (maxx - minx) * (maxy - miny)

        if area < best_area:
            rect_rot = Polygon([
                (minx, miny),
                (maxx, miny),
                (maxx, maxy),
                (minx, maxy),
                (minx, miny),
            ])

            # Rotate rectangle back to original coordinates
            rect = affinity.rotate(rect_rot, angle, origin="centroid", use_radians=False)

            best_area = area
            best_rect = rect
            best_angle = angle

    xy = np.asarray(best_rect.exterior.coords)

    return xy, best_angle, best_rect


def clip_elements_to_polygon(elem_gdf, clip_poly, clip_elem=True):
    """
    Select element polygons overlapping a model domain polygon and
    optionally clip them to the domain boundary.

    Parameters
    ----------
    elem_gdf : GeoDataFrame
        Input element polygons.
    clip_poly : shapely.geometry.Polygon
        Polygon defining the local model domain.
    clip_elem : bool, default=True
        If True, clip intersecting elements to the domain.
        If False, return the original geometries of all
        intersecting elements.
    """

    # Determine relationship to model domain
    inside = elem_gdf.within(clip_poly)
    intersect = elem_gdf.intersects(clip_poly)

    # Keep only intersecting elements
    work_gdf = elem_gdf.loc[intersect].copy()

    # Add flags
    work_gdf["is_internal"] = inside.loc[intersect].values
    work_gdf["is_boundary_clipped"] = (~inside.loc[intersect]).values

    if not clip_elem:
        return work_gdf.reset_index(drop=True)

    # Clip geometries
    rect_gdf = gpd.GeoDataFrame(
        geometry=[clip_poly],
        crs=elem_gdf.crs
    )

    clipped_gdf = gpd.clip(work_gdf, rect_gdf)

    return clipped_gdf.reset_index(drop=True)

def isolate_submodel_mesh(elem_df, node_df, submodel_elem,
    elem_dict=None, node_dict=None, submodel_dict=None):
    """
    Extract nodes associated with a submodel element subset and create
    local node/element indexing.

    Returns
    -------
    sub_node_df : DataFrame
        Node subset with fields:
        IDnew, original node ID, and all original node fields.

    sub_elem_df : DataFrame
        Element subset with fields:
        IEnew, original IE, all original element fields,
        NDinew columns, and original NDi columns.
    """

    # Defaults
    elem_dict = elem_dict or {}
    node_dict = node_dict or {}
    submodel_dict = submodel_dict or {}

    elem_id_col = elem_dict.get("elem_id", "IE")
    elem_node_cols = elem_dict.get("node_cols", ["ND1", "ND2", "ND3", "ND4"])

    node_id_col = node_dict.get("node_id", "ID")

    sub_elem_id_col = submodel_dict.get("elem_id", "ElementID")

    # Keep only node columns that actually exist
    elem_node_cols = [c for c in elem_node_cols if c in elem_df.columns]

    # Merge submodel fields with full element definition
    sub_elem_df = submodel_elem.merge(
        elem_df,
        left_on=sub_elem_id_col,
        right_on=elem_id_col,
        how="inner",
        suffixes=("_submodel", "")
    ).copy()

    # Add new element index
    sub_elem_df.insert(0, "IEnew", np.arange(1, len(sub_elem_df) + 1))

    # Collect used node IDs
    node_ids = sub_elem_df[elem_node_cols].to_numpy().ravel()
    node_ids = pd.Series(node_ids)
    node_ids = node_ids[(node_ids.notna()) & (node_ids != 0)].unique()

    # Keep all original node fields
    sub_node_df = node_df.loc[
        node_df[node_id_col].isin(node_ids)
    ].copy()

    sub_node_df.insert(0, "IDnew", np.arange(1, len(sub_node_df) + 1))

    # Original node ID -> local node ID
    node_id_to_new = dict(zip(sub_node_df[node_id_col], sub_node_df["IDnew"]))

    # Add local node columns
    for nd_col in elem_node_cols:
        sub_elem_df[f"{nd_col}new"] = (
            sub_elem_df[nd_col]
            .map(node_id_to_new)
            .fillna(0)
            .astype(int)
        )

    return sub_node_df.reset_index(drop=True), sub_elem_df.reset_index(drop=True)


def extract_clipped_mesh_outline(
    elem_clipped,
    elem_id_col="ElementID",
    attrs=None,
    round_decimals=8,
):
    """
    Extract ordered outline points and boundary segments from a clipped mesh.

    Parameters
    ----------
    elem_clipped : GeoDataFrame
        Clipped element polygons. Should contain ElementID and geometry.
    elem_id_col : str, default="ElementID"
        Element ID field.
    attrs : list[str] or None
        Additional element attributes to carry to the segment dataframe.
    round_decimals : int, default=8
        Decimal precision used to identify matching shared edges.

    Returns
    -------
    outline_points_df : DataFrame
        Ordered boundary points with fields:
        ring_id, point_id, x, y

    outline_segments_df : DataFrame
        Ordered boundary segments with fields:
        ring_id, segment_id,
        p0, p1,
        x0, y0, x1, y1,
        ElementID,
        plus requested attrs.
    """

    if attrs is None:
        attrs = []

    attrs = [a for a in attrs if a in elem_clipped.columns]

    def key_xy(x, y):
        return (round(float(x), round_decimals),
                round(float(y), round_decimals))

    def polygon_rings(geom):
        """Return exterior and interior rings from Polygon/MultiPolygon."""
        if geom.geom_type == "Polygon":
            yield list(geom.exterior.coords)
            for interior in geom.interiors:
                yield list(interior.coords)

        elif geom.geom_type == "MultiPolygon":
            for poly in geom.geoms:
                yield list(poly.exterior.coords)
                for interior in poly.interiors:
                    yield list(interior.coords)

    # ------------------------------------------------------------------
    # 1. Collect all element edges.
    #    Shared internal edges appear twice.
    #    External boundary edges appear once.
    # ------------------------------------------------------------------
    edge_records = defaultdict(list)

    for _, row in elem_clipped.iterrows():
        geom = row.geometry

        if geom is None or geom.is_empty:
            continue

        for coords in polygon_rings(geom):
            for i in range(len(coords) - 1):
                x0, y0 = coords[i][:2]
                x1, y1 = coords[i + 1][:2]

                k0 = key_xy(x0, y0)
                k1 = key_xy(x1, y1)

                if k0 == k1:
                    continue

                edge_key = tuple(sorted([k0, k1]))

                rec = {
                    "k0": k0,
                    "k1": k1,
                    "x0": float(x0),
                    "y0": float(y0),
                    "x1": float(x1),
                    "y1": float(y1),
                    elem_id_col: row[elem_id_col],
                }

                for a in attrs:
                    rec[a] = row[a]

                edge_records[edge_key].append(rec)

    # Keep only edges that occur once
    boundary_edges = [
        records[0]
        for records in edge_records.values()
        if len(records) == 1
    ]

    # ------------------------------------------------------------------
    # 2. Build node-edge adjacency for ordering.
    # ------------------------------------------------------------------
    adjacency = defaultdict(list)

    for i, e in enumerate(boundary_edges):
        adjacency[e["k0"]].append(i)
        adjacency[e["k1"]].append(i)

    unused = set(range(len(boundary_edges)))

    point_rows = []
    segment_rows = []

    ring_id = 0

    while unused:
        ring_id += 1

        # start from any unused edge
        start_edge_idx = next(iter(unused))
        e0 = boundary_edges[start_edge_idx]

        current = e0["k0"]
        start_node = current

        ordered_points = [current]
        ordered_edges = []

        while True:
            candidate_edges = [
                ei for ei in adjacency[current]
                if ei in unused
            ]

            if not candidate_edges:
                break

            ei = candidate_edges[0]
            edge = boundary_edges[ei]
            unused.remove(ei)

            k0 = edge["k0"]
            k1 = edge["k1"]

            if current == k0:
                next_node = k1
                oriented = edge.copy()
            else:
                next_node = k0
                oriented = edge.copy()
                oriented["x0"], oriented["x1"] = oriented["x1"], oriented["x0"]
                oriented["y0"], oriented["y1"] = oriented["y1"], oriented["y0"]
                oriented["k0"], oriented["k1"] = oriented["k1"], oriented["k0"]

            ordered_edges.append(oriented)
            ordered_points.append(next_node)
            current = next_node

            if current == start_node:
                break

        # Write points
        point_offset = len(point_rows)

        for ip, k in enumerate(ordered_points):
            point_rows.append({
                "ring_id": ring_id,
                "point_id": ip,
                "x": k[0],
                "y": k[1],
            })

        # Write segments
        for iseg, edge in enumerate(ordered_edges):
            row = {
                "ring_id": ring_id,
                "segment_id": iseg,
                "p0": point_offset + iseg,
                "p1": point_offset + iseg + 1,
                "x0": edge["x0"],
                "y0": edge["y0"],
                "x1": edge["x1"],
                "y1": edge["y1"],
                elem_id_col: edge[elem_id_col],
            }

            for a in attrs:
                row[a] = edge[a]

            segment_rows.append(row)

    outline_points_df = pd.DataFrame(point_rows)
    outline_segments_df = pd.DataFrame(segment_rows)

    return outline_points_df, outline_segments_df