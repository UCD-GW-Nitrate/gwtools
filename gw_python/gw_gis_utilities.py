import pandas as pd
import geopandas as gpd
import numpy as np
from scipy.spatial import Delaunay
from pathlib import Path

from shapely.geometry import (
    Point, LineString, Polygon, MultiPolygon,
    MultiLineString, GeometryCollection
)
from shapely.ops import unary_union, linemerge
import plotly.graph_objects as go

def assign_fields_from_polygons_to_points(points: pd.DataFrame, polygons: gpd.GeoDataFrame,
                                          fields_to_copy: list, point_crs: str = "EPSG:3310",
                                          unique: bool = False) -> gpd.GeoDataFrame:
    """
        Select points from A that fall inside polygons from B,
        and append specified fields from B. Returns a GeoDataFrame.

        Parameters
        ----------
        A : pd.DataFrame
            Must have columns ['X', 'Y'].
        B : gpd.GeoDataFrame
            Polygon layer.
        C : list
            List of column names from B to append.
        crs : str, optional
            Coordinate reference system for the points in A. Default "EPSG:4326".
        unique : bool, optional
            If True, keep only the first polygon match per point.
            If False (default), return all matches.

        Returns
        -------
        gpd.GeoDataFrame
            Subset of A inside polygons of B, with appended fields.
        """

    # Convert points to GeoDataFrame
    gdf_points = gpd.GeoDataFrame(
        points.copy(),
        geometry=gpd.points_from_xy(points['X'], points['Y']),
        crs=point_crs
    )

    # Ensure CRS matches
    if gdf_points.crs != polygons.crs:
        gdf_points = gdf_points.to_crs(polygons.crs)

    # Spatial join (points within polygons)
    joined = gpd.sjoin(gdf_points, polygons[fields_to_copy + ["geometry"]], how="inner", predicate="within")
    # Drop unnecessary column
    joined = joined.drop(columns=["index_right"])

    # If unique requested → keep first match per point
    if unique:
        joined = joined.drop_duplicates(subset=gdf_points.columns.tolist())

    return joined.reset_index(drop=True)

def triangulate_inside_outline(outline_gdf, node_gdf):
    """
    Keep point nodes inside outline_gdf, triangulate them, then remove
    triangles whose barycenters fall outside outline_gdf.

    Returns
    -------
    node_gdf_new : geopandas.GeoDataFrame
        Filtered node GeoDataFrame with original fields preserved.

    tri_ids_new : np.ndarray
        Triangle connectivity using the new node_gdf_new row order,
        shape (Ntri, 3).
    """

    outline_geom = outline_gdf.geometry.union_all()

    # Keep points inside or on outline
    keep_nodes = (
        node_gdf.geometry.within(outline_geom) |
        node_gdf.geometry.touches(outline_geom)
    )

    node_gdf_new = node_gdf.loc[keep_nodes].copy().reset_index(drop=True)

    if len(node_gdf_new) < 3:
        return node_gdf_new, np.empty((0, 3), dtype=int)

    # Extract X/Y from point geometry, including POINT Z
    xy_nodes = np.column_stack([
        node_gdf_new.geometry.x.to_numpy(),
        node_gdf_new.geometry.y.to_numpy()
    ])

    # Triangulate retained nodes
    tri = Delaunay(xy_nodes)
    tri_ids_all = tri.simplices.copy()

    # Remove triangles with barycenters outside outline
    barycenters = xy_nodes[tri_ids_all].mean(axis=1)

    bary_points = gpd.GeoSeries(
        gpd.points_from_xy(barycenters[:, 0], barycenters[:, 1]),
        crs=node_gdf.crs
    )

    keep_tri = (
        bary_points.within(outline_geom) |
        bary_points.touches(outline_geom)
    )

    tri_ids_new = tri_ids_all[np.asarray(keep_tri)]

    return node_gdf_new, tri_ids_new


def write_gpkg_layer(gdf: gpd.GeoDataFrame,
                     gpkg_name: str,
                     layer_name: str) -> None:
    """
    Write a GeoDataFrame to a GeoPackage layer.

    If the GeoPackage does not exist, it is created.
    If the layer already exists, it is overwritten.
    Other layers in the GeoPackage are preserved.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame to write.

    gpkg_name : str
        Path to the GeoPackage file (*.gpkg).

    layer_name : str
        Name of the layer within the GeoPackage.

    Examples
    --------
    write_gpkg_layer(outline_gdf,
                     "model_layers.gpkg",
                     "outline")

    write_gpkg_layer(outline_gdf_buff,
                     "model_layers.gpkg",
                     "outline_buff")
    """
    gpkg_name = str(Path(gpkg_name))

    gdf.to_file(
        gpkg_name,
        layer=layer_name,
        driver="GPKG",
        mode="w"
    )

def project_boundary_parts_to_npsat(
        parts,
        msh_outline,
        max_distance=None,
        n_candidates=20,
        endpoint_tol=1e-9,
):
    """
    Project C2VSim boundary-part nodes onto NPSAT boundary segments.

    Projection criterion:
        1. Transform parts to the CRS of msh_outline.
        2. Find nearby NPSAT boundary segments.
        3. Accept only segments where the perpendicular
           projection falls within the segment.
        4. Choose the candidate with minimum perpendicular distance.

    Parameters
    ----------
    parts : list of dict
        Output of get_stable_boundary_parts().
        Each dictionary contains 'nodes', 'lines', etc.

    msh_outline : GeoDataFrame
        NPSAT mesh outline (EPSG:3310).
        Can contain polygons or lines.

    max_distance : float, optional
        Maximum acceptable projection distance (m).
        If exceeded, projected coordinates are NaN.

    n_candidates : int
        Initial number of nearby segments examined.
        The search expands if no valid projection is found.

    endpoint_tol : float
        Numerical tolerance for projected position along a segment,
        expressed as a fraction of segment length.

    Returns
    -------
    projected_parts : list of dict
        Copies of the input parts with additional node fields:

        proj_x       : projected X coordinate (EPSG:3310)
        proj_y       : projected Y coordinate (EPSG:3310)
        proj_dist    : projection distance (m)
        proj_seg_id  : NPSAT boundary segment index

    Notes
    -----
    Original node geometry is transformed to EPSG:3310,
    but not replaced with the projected geometry.

    The original 'X', 'Y' columns, if present, are retained.
    """

    if msh_outline.crs is None:
        raise ValueError("msh_outline must have a defined CRS.")

    if n_candidates < 1:
        raise ValueError("n_candidates must be >= 1.")

    target_crs = msh_outline.crs

    # -------------------------------------------------------
    # 1. Extract NPSAT outline segments
    # -------------------------------------------------------
    segments = []

    def add_segments(geom):

        if geom is None or geom.is_empty:
            return

        if isinstance(geom, Polygon):
            add_segments(geom.boundary)

        elif isinstance(geom, MultiPolygon):
            for g in geom.geoms:
                add_segments(g.boundary)

        elif isinstance(geom, LineString):
            coords = list(geom.coords)

            for p0, p1 in zip(coords[:-1], coords[1:]):
                if p0 != p1:
                    segments.append(LineString([p0, p1]))

        elif isinstance(geom, (MultiLineString, GeometryCollection)):
            for g in geom.geoms:
                add_segments(g)

    # Dissolve polygons to avoid internal mesh edges
    polygon_geoms = []
    line_geoms = []

    for geom in msh_outline.geometry:
        if geom is None or geom.is_empty:
            continue

        if isinstance(geom, (Polygon, MultiPolygon)):
            polygon_geoms.append(geom)
        else:
            line_geoms.append(geom)

    if polygon_geoms:
        add_segments(unary_union(polygon_geoms).boundary)

    for geom in line_geoms:
        add_segments(geom)

    if not segments:
        raise ValueError("No valid NPSAT boundary segments found.")

    seg_gdf = gpd.GeoDataFrame(
        {"seg_id": np.arange(len(segments))},
        geometry=segments,
        crs=target_crs,
    )

    # Segment endpoint arrays
    xy0 = np.array([s.coords[0] for s in segments], dtype=float)[:, :2]
    xy1 = np.array([s.coords[-1] for s in segments], dtype=float)[:, :2]

    vec = xy1 - xy0
    length2 = np.sum(vec ** 2, axis=1)

    sindex = seg_gdf.sindex
    n_segments = len(segments)

    # -------------------------------------------------------
    # 2. Project individual point
    # -------------------------------------------------------
    def project_point(x, y):

        p = Point(x, y)

        # Nearest candidates using spatial index
        k = min(n_candidates, n_segments)

        while True:

            # Query nearby segments through expanding search radius
            nearest_idx = list(
                sindex.nearest(
                    p,
                    return_all=True
                )[1]
            )

            # Start with nearest segment(s)
            candidate_ids = set(nearest_idx)

            if len(candidate_ids) < k:
                # Include additional nearby segments by
                # expanding a spatial search envelope.
                radius = max(
                    min(segments[j].distance(p) for j in candidate_ids),
                    1.0
                )

                while len(candidate_ids) < k:
                    radius *= 2.0

                    candidate_ids.update(
                        sindex.query(p.buffer(radius))
                    )

                    if len(candidate_ids) == n_segments:
                        break

            ids = np.fromiter(candidate_ids, dtype=int)

            # Vectorized perpendicular projection
            w = np.array([x, y]) - xy0[ids]

            t = np.sum(w * vec[ids], axis=1) / length2[ids]

            # Projection must fall within the segment
            valid = (
                    (t >= -endpoint_tol)
                    & (t <= 1.0 + endpoint_tol)
            )

            if valid.any():

                valid_ids = ids[valid]
                valid_t = np.clip(t[valid], 0.0, 1.0)

                proj = (
                        xy0[valid_ids]
                        + valid_t[:, None] * vec[valid_ids]
                )

                distances = np.linalg.norm(
                    proj - np.array([x, y]),
                    axis=1
                )

                j = np.argmin(distances)

                if (
                        max_distance is not None
                        and distances[j] > max_distance
                ):
                    return np.nan, np.nan, np.nan, -1

                return (
                    proj[j, 0],
                    proj[j, 1],
                    distances[j],
                    int(valid_ids[j]),
                )

            if k == n_segments:
                break

            k = min(2 * k, n_segments)

        return np.nan, np.nan, np.nan, -1

    # -------------------------------------------------------
    # 3. Process all boundary parts
    # -------------------------------------------------------
    projected_parts = []

    for part in parts:
        new_part = part.copy()

        nodes = part["nodes"].to_crs(target_crs).copy()

        coords = np.column_stack([
            nodes.geometry.x.to_numpy(),
            nodes.geometry.y.to_numpy()
        ])

        results = np.array(
            [project_point(x, y) for x, y in coords]
        )

        nodes["proj_x"] = results[:, 0]
        nodes["proj_y"] = results[:, 1]
        nodes["proj_dist"] = results[:, 2]
        nodes["proj_seg_id"] = results[:, 3].astype(int)

        new_part["nodes"] = nodes

        # Also transform original boundary lines
        new_part["lines"] = part["lines"].to_crs(target_crs)

        projected_parts.append(new_part)

    return projected_parts


def plot_boundary_projection(
    projected_parts,
    msh_outline,
    original_color="green",
    projection_line_color="orange",
    marker_size=6,
    fig=None,
):
    """
    Visualize C2VSim boundary nodes projected onto the NPSAT outline.

    Parameters
    ----------
    projected_parts : list of dict
        Output of project_boundary_parts_to_npsat().
        Each part contains a 'nodes' GeoDataFrame with:
        proj_x, proj_y, proj_dist, ID.

    msh_outline : GeoDataFrame
        NPSAT mesh outline. Automatically transformed to
        the CRS of projected_parts if necessary.

    original_color : str
        Color of original C2VSim boundary nodes.

    projection_line_color : str
        Color of projection lines.

    marker_size : float
        Marker size for original and projected nodes.

    fig : go.Figure, optional
        Existing figure.

    Returns
    -------
    fig : go.Figure
    """

    if fig is None:
        fig = go.Figure()

    if not projected_parts:
        raise ValueError("projected_parts is empty.")

    target_crs = projected_parts[0]["nodes"].crs

    if target_crs is None or msh_outline.crs is None:
        raise ValueError("Both datasets must have a defined CRS.")

    outline = msh_outline.to_crs(target_crs)

    # -------------------------------------------------------
    # 1. NPSAT outline: one trace
    # -------------------------------------------------------
    outline_x, outline_y = [], []

    def add_geometry(geom):

        if geom is None or geom.is_empty:
            return

        if isinstance(geom, Polygon):
            add_geometry(geom.boundary)

        elif isinstance(geom, MultiPolygon):
            for g in geom.geoms:
                add_geometry(g.boundary)

        elif isinstance(geom, LineString):
            x, y = geom.xy

            outline_x.extend(x)
            outline_x.append(None)

            outline_y.extend(y)
            outline_y.append(None)

        elif isinstance(geom, (MultiLineString, GeometryCollection)):
            for g in geom.geoms:
                add_geometry(g)

    for geom in outline.geometry:
        add_geometry(geom)

    fig.add_trace(go.Scatter(
        x=outline_x,
        y=outline_y,
        mode="lines",
        line=dict(color="gray", width=1.5),
        name="NPSAT boundary",
        hoverinfo="skip",
        showlegend=True,
    ))

    # -------------------------------------------------------
    # 2. Combine nodes from all parts
    # -------------------------------------------------------
    frames = []

    for i, part in enumerate(projected_parts):

        nodes = part["nodes"].copy()

        nodes["part_id"] = i + 1
        nodes["part_length_km"] = part["length_km"]

        frames.append(nodes)

    nodes = gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True),
        geometry="geometry",
        crs=target_crs,
    )

    valid = (
        np.isfinite(nodes["proj_x"].to_numpy(dtype=float))
        & np.isfinite(nodes["proj_y"].to_numpy(dtype=float))
        & np.isfinite(nodes["proj_dist"].to_numpy(dtype=float))
    )

    nodes = nodes.loc[valid].copy()

    if nodes.empty:
        raise ValueError("No successfully projected nodes found.")

    x0 = nodes.geometry.x.to_numpy(dtype=float)
    y0 = nodes.geometry.y.to_numpy(dtype=float)

    x1 = nodes["proj_x"].to_numpy(dtype=float)
    y1 = nodes["proj_y"].to_numpy(dtype=float)

    dist = nodes["proj_dist"].to_numpy(dtype=float)
    node_id = nodes["ID"].to_numpy()
    part_id = nodes["part_id"].to_numpy()

    # -------------------------------------------------------
    # 3. Projection lines: one trace
    # -------------------------------------------------------
    # Each line has:
    # original point -> projected point -> None
    n = len(nodes)

    line_x = np.empty(3 * n, dtype=object)
    line_y = np.empty(3 * n, dtype=object)

    line_x[0::3] = x0
    line_x[1::3] = x1
    line_x[2::3] = None

    line_y[0::3] = y0
    line_y[1::3] = y1
    line_y[2::3] = None

    # Repeat metadata at both ends of each projection line
    line_customdata = np.empty((3 * n, 3), dtype=object)

    line_customdata[0::3, 0] = part_id
    line_customdata[1::3, 0] = part_id

    line_customdata[0::3, 1] = node_id
    line_customdata[1::3, 1] = node_id

    line_customdata[0::3, 2] = dist
    line_customdata[1::3, 2] = dist

    line_customdata[2::3, :] = None

    fig.add_trace(go.Scatter(
        x=line_x.tolist(),
        y=line_y.tolist(),
        mode="lines",
        line=dict(
            color=projection_line_color,
            width=1,
        ),
        customdata=line_customdata,
        hovertemplate=(
            "<b>Projection line</b><br>"
            "Part ID: %{customdata[0]}<br>"
            "Node ID: %{customdata[1]}<br>"
            "Distance: %{customdata[2]:.2f} m"
            "<extra></extra>"
        ),
        name="Projection lines",
        showlegend=False,
    ))

    # -------------------------------------------------------
    # 4. Original nodes: one trace
    # -------------------------------------------------------
    point_customdata = np.column_stack([
        part_id,
        node_id,
        dist,
    ])

    fig.add_trace(go.Scatter(
        x=x0,
        y=y0,
        mode="markers",
        marker=dict(
            size=marker_size,
            color=original_color,
            symbol="circle",
        ),
        customdata=point_customdata,
        hovertemplate=(
            "<b>Original C2VSim node</b><br>"
            "Part ID: %{customdata[0]}<br>"
            "Node ID: %{customdata[1]}<br>"
            "Projection distance: %{customdata[2]:.2f} m"
            "<extra></extra>"
        ),
        name="Original nodes",
        showlegend=False,
    ))

    # -------------------------------------------------------
    # 5. Projected nodes: one trace, colored by distance
    # -------------------------------------------------------
    fig.add_trace(go.Scatter(
        x=x1,
        y=y1,
        mode="markers",
        marker=dict(
            size=marker_size,
            color=dist,
            colorscale="Viridis",
            cmin=0,
            cmax=max(float(np.max(dist)), 1e-9),
            showscale=True,
            colorbar=dict(title="Projection (m)"),
            symbol="circle",
            line=dict(width=0.5, color="black"),
        ),
        customdata=point_customdata,
        hovertemplate=(
            "<b>Projected NPSAT node</b><br>"
            "Part ID: %{customdata[0]}<br>"
            "Original node ID: %{customdata[1]}<br>"
            "Projection distance: %{customdata[2]:.2f} m"
            "<extra></extra>"
        ),
        name="Projected nodes",
        showlegend=False,
    ))

    # -------------------------------------------------------
    # 6. Figure layout
    # -------------------------------------------------------
    fig.update_layout(
        title=(
            "C2VSim → NPSAT boundary projection"
            f"<br><sup>{n} projected nodes | "
            f"{len(projected_parts)} boundary parts | "
            f"Mean distance: {np.mean(dist):.1f} m | "
            f"Max distance: {np.max(dist):.1f} m</sup>"
        ),
        xaxis_title="Easting (m)",
        yaxis_title="Northing (m)",
        #template="plotly_white",
        width=1100,
        height=850,
        hovermode="closest",
    )

    fig.update_yaxes(
        scaleanchor="x",
        scaleratio=1,
    )

    return fig


def reconstruct_npsat_boundary_parts(
    bc_parts,
    bc_heads,
    bc_elev,
    msh_outline,
    snap_tolerance=10.0,
    max_projection_distance=20.0,
):
    """
    Reconstruct boundary-condition parts along the actual NPSAT
    mesh outline, inserting intermediate mesh vertices and
    interpolating heads and elevations.

    Parameters
    ----------
    bc_parts : list of dict
        Boundary parts containing 'nodes' GeoDataFrames.
        Required node fields:
            row_id, proj_x, proj_y

        Optional fields:
            ID, original_part_id, subpart_id

    bc_heads : ndarray, shape (N, n_layers)
        Mean groundwater heads.

    bc_elev : ndarray, shape (N, n_layers - 1)
        Interlayer elevations.

    msh_outline : GeoDataFrame
        NPSAT mesh outline, typically EPSG:3310.

    snap_tolerance : float
        Distance (m) within which projected nodes are snapped
        to actual NPSAT mesh vertices.

    max_projection_distance : float
        Maximum allowed distance between a projected point
        and the NPSAT outline.

    Returns
    -------
    final_parts : list of dict
        Reconstructed boundary parts with:
            nodes : GeoDataFrame
            lines : GeoDataFrame
            length
            length_km
            n_nodes
            n_segments

    final_heads : ndarray
        Interpolated heads, indexed by nodes.row_id.

    final_elev : ndarray
        Interpolated interlayer elevations, indexed by nodes.row_id.

    Notes
    -----
    - The NPSAT outline must represent one connected boundary loop.
    - Input nodes must be ordered along their original boundary parts.
    - New nodes receive ID = NaN.
    - The returned geometry is in msh_outline.crs.
    - No interpolation is performed across gaps between parts.
    """

    # -------------------------------------------------------
    # 0. Validate inputs
    # -------------------------------------------------------
    bc_heads = np.asarray(bc_heads, dtype=float)
    bc_elev = np.asarray(bc_elev, dtype=float)

    if bc_heads.ndim != 2 or bc_elev.ndim != 2:
        raise ValueError("bc_heads and bc_elev must be 2D arrays.")

    if bc_heads.shape[0] != bc_elev.shape[0]:
        raise ValueError("Head and elevation row counts differ.")

    if msh_outline.crs is None:
        raise ValueError("msh_outline must have a defined CRS.")

    if snap_tolerance < 0 or max_projection_distance < 0:
        raise ValueError("Tolerances must be nonnegative.")

    # -------------------------------------------------------
    # 1. Extract a single ordered NPSAT boundary ring
    # -------------------------------------------------------
    geometries = [
        g for g in msh_outline.geometry
        if g is not None and not g.is_empty
    ]

    if not geometries:
        raise ValueError("NPSAT outline is empty.")

    if all(isinstance(g, (Polygon, MultiPolygon)) for g in geometries):
        boundary = unary_union(geometries).boundary
    else:
        boundary = unary_union([
            g.boundary if isinstance(g, (Polygon, MultiPolygon))
            else g
            for g in geometries
        ])

    if isinstance(boundary, MultiLineString):
        boundary = linemerge(boundary)

    if not isinstance(boundary, LineString):
        raise ValueError(
            "Expected one connected NPSAT boundary LineString. "
            "Check for multiple loops or disconnected geometry."
        )

    coords = np.asarray(boundary.coords, dtype=float)[:, :2]

    if not np.allclose(coords[0], coords[-1], atol=1e-8, rtol=0):
        raise ValueError("NPSAT boundary must be a closed loop.")

    # Remove duplicate consecutive vertices
    step = np.linalg.norm(np.diff(coords, axis=0), axis=1)

    keep = np.r_[True, step > 1e-9]
    coords = coords[keep]

    if not np.allclose(coords[0], coords[-1], atol=1e-8, rtol=0):
        coords = np.vstack([coords, coords[0]])

    # Boundary cumulative distances
    segment_lengths = np.linalg.norm(
        np.diff(coords, axis=0), axis=1
    )

    mesh_s = np.r_[0.0, np.cumsum(segment_lengths)]
    perimeter = mesh_s[-1]

    boundary_line = LineString(coords)

    # Unique mesh vertices (exclude repeated closing vertex)
    mesh_vertices = coords[:-1]
    mesh_vertex_s = mesh_s[:-1]

    # -------------------------------------------------------
    # 2. Helper: locate point on mesh outline
    # -------------------------------------------------------
    def locate_on_boundary(x, y):

        p = Point(float(x), float(y))

        s = boundary_line.project(p)
        q = boundary_line.interpolate(s)

        projected_xy = np.array([q.x, q.y])

        residual = np.linalg.norm(
            projected_xy - np.array([x, y])
        )

        if residual > max_projection_distance:
            raise ValueError(
                f"Point ({x:.2f}, {y:.2f}) is "
                f"{residual:.2f} m from NPSAT boundary."
            )

        # Find closest mesh vertex
        d = np.linalg.norm(
            mesh_vertices - projected_xy,
            axis=1
        )

        nearest_vertex = int(np.argmin(d))

        if d[nearest_vertex] <= snap_tolerance:

            return (
                mesh_vertex_s[nearest_vertex],
                mesh_vertices[nearest_vertex].copy(),
                nearest_vertex,
                "snapped",
            )

        return (
            float(s),
            projected_xy,
            -1,
            "original",
        )

    # -------------------------------------------------------
    # 3. Helper: unwrap boundary distances
    # -------------------------------------------------------
    def unwrap_distances(s_values):

        """
        Convert periodic boundary coordinates into a
        continuous ordered sequence.

        Chooses forward or backward orientation based on
        the shorter total path along the boundary.
        """

        s_values = np.asarray(s_values, dtype=float)

        if len(s_values) < 2:
            return s_values.copy()

        forward = np.mod(
            np.diff(s_values), perimeter
        )

        backward = np.mod(
            -np.diff(s_values), perimeter
        )

        if forward.sum() <= backward.sum():
            direction = 1
            increments = forward
        else:
            direction = -1
            increments = backward

        unwrapped = np.r_[
            0.0,
            np.cumsum(increments)
        ]

        return unwrapped, direction

    # -------------------------------------------------------
    # 4. Reconstruct each part
    # -------------------------------------------------------
    final_parts = []
    all_heads = []
    all_elev = []

    row_offset = 0

    for part_idx, part in enumerate(bc_parts):

        nodes = part["nodes"].copy().reset_index(drop=True)

        if len(nodes) < 2:
            continue

        if nodes.crs != msh_outline.crs:
            nodes = nodes.to_crs(msh_outline.crs)

        old_row_ids = nodes["row_id"].to_numpy(dtype=int)

        part_heads = bc_heads[old_row_ids]
        part_elev = bc_elev[old_row_ids]

        # ---------------------------------------------------
        # Locate and snap original projected nodes
        # ---------------------------------------------------
        located = [
            locate_on_boundary(x, y)
            for x, y in nodes[
                ["proj_x", "proj_y"]
            ].to_numpy(dtype=float)
        ]

        s_values = np.array([v[0] for v in located])
        xy_values = np.array([v[1] for v in located])
        vertex_ids = np.array([v[2] for v in located])
        sources = [v[3] for v in located]

        # ---------------------------------------------------
        # Unwrap cumulative boundary distances
        # ---------------------------------------------------
        part_s, direction = unwrap_distances(s_values)

        # ---------------------------------------------------
        # Remove duplicate consecutive anchor positions
        # ---------------------------------------------------
        keep_anchor = np.r_[
            True,
            np.diff(part_s) > 1e-8
        ]

        part_s = part_s[keep_anchor]
        xy_values = xy_values[keep_anchor]
        vertex_ids = vertex_ids[keep_anchor]

        part_heads = part_heads[keep_anchor]
        part_elev = part_elev[keep_anchor]

        nodes = nodes.loc[keep_anchor].reset_index(drop=True)

        sources = [
            src for src, keep in zip(sources, keep_anchor)
            if keep
        ]

        if len(part_s) < 2:
            continue

        # ---------------------------------------------------
        # Build reconstructed nodes
        # ---------------------------------------------------
        records = []
        new_heads = []
        new_elev = []

        def append_node(
            x, y, s_local,
            node_id, mesh_vertex_id,
            source, head, elev
        ):

            records.append({
                "ID": node_id,
                "X": x,
                "Y": y,
                "proj_x": x,
                "proj_y": y,
                "part_distance": s_local,
                "mesh_vertex_id": mesh_vertex_id,
                "node_source": source,
                "geometry": Point(x, y),
            })

            new_heads.append(head)
            new_elev.append(elev)

        # First anchor
        append_node(
            xy_values[0, 0],
            xy_values[0, 1],
            part_s[0],
            nodes.iloc[0].get("ID", np.nan),
            vertex_ids[0],
            sources[0],
            part_heads[0],
            part_elev[0],
        )

        # ---------------------------------------------------
        # Process consecutive anchors
        # ---------------------------------------------------
        for i in range(len(part_s) - 1):

            s0 = part_s[i]
            s1 = part_s[i + 1]

            if s1 <= s0:
                continue

            # Boundary coordinate of first anchor
            mesh_start = s_values[keep_anchor][i]

            # Determine intermediate mesh vertices
            intermediate = []

            for j, vertex_s in enumerate(mesh_vertex_s):

                if direction == 1:
                    delta = (vertex_s - mesh_start) % perimeter
                else:
                    delta = (mesh_start - vertex_s) % perimeter

                if 1e-8 < delta < (s1 - s0) - 1e-8:
                    intermediate.append((delta, j))

            intermediate.sort(key=lambda item: item[0])

            # Insert mesh vertices with interpolated data
            for delta, j in intermediate:

                weight = delta / (s1 - s0)

                head = (
                    (1 - weight) * part_heads[i]
                    + weight * part_heads[i + 1]
                )

                elev = (
                    (1 - weight) * part_elev[i]
                    + weight * part_elev[i + 1]
                )

                append_node(
                    mesh_vertices[j, 0],
                    mesh_vertices[j, 1],
                    s0 + delta,
                    np.nan,
                    j,
                    "interpolated",
                    head,
                    elev,
                )

            # Append next anchor
            append_node(
                xy_values[i + 1, 0],
                xy_values[i + 1, 1],
                s1,
                nodes.iloc[i + 1].get("ID", np.nan),
                vertex_ids[i + 1],
                sources[i + 1],
                part_heads[i + 1],
                part_elev[i + 1],
            )

        # ---------------------------------------------------
        # Build final GeoDataFrames
        # ---------------------------------------------------
        final_nodes = gpd.GeoDataFrame(
            records,
            geometry="geometry",
            crs=msh_outline.crs,
        )

        final_heads_part = np.asarray(new_heads)
        final_elev_part = np.asarray(new_elev)

        # Remove any zero-length consecutive segments
        xy = final_nodes[
            ["X", "Y"]
        ].to_numpy(dtype=float)

        lengths = np.linalg.norm(
            np.diff(xy, axis=0),
            axis=1
        )

        keep_nodes = np.r_[
            True,
            lengths > 1e-8
        ]

        final_nodes = (
            final_nodes.loc[keep_nodes]
            .reset_index(drop=True)
        )

        final_heads_part = final_heads_part[keep_nodes]
        final_elev_part = final_elev_part[keep_nodes]

        if len(final_nodes) < 2:
            continue

        xy = final_nodes[["X", "Y"]].to_numpy()

        lengths = np.linalg.norm(
            np.diff(xy, axis=0), axis=1
        )

        final_nodes["part_distance"] = np.r_[
            0.0,
            np.cumsum(lengths)
        ]

        n_final = len(final_nodes)

        final_nodes["row_id"] = np.arange(
            row_offset,
            row_offset + n_final,
            dtype=int
        )

        # Reconstructed boundary segments
        line_geoms = [
            LineString([xy[j], xy[j + 1]])
            for j in range(n_final - 1)
        ]

        final_lines = gpd.GeoDataFrame(
            {
                "length": lengths,
                "part_distance": np.r_[
                    0.0,
                    np.cumsum(lengths[:-1])
                ],
            },
            geometry=line_geoms,
            crs=msh_outline.crs,
        )

        total_length = float(lengths.sum())

        new_part = part.copy()

        new_part.update({
            "nodes": final_nodes,
            "lines": final_lines,
            "length": total_length,
            "length_km": total_length / 1000.0,
            "n_nodes": n_final,
            "n_segments": len(final_lines),
            "original_part_id": part.get(
                "original_part_id", part_idx + 1
            ),
            "subpart_id": part.get("subpart_id", 1),
        })

        final_parts.append(new_part)

        all_heads.append(final_heads_part)
        all_elev.append(final_elev_part)

        row_offset += n_final

    # -------------------------------------------------------
    # 5. Assemble global arrays
    # -------------------------------------------------------
    final_heads = (
        np.vstack(all_heads)
        if all_heads
        else np.empty((0, bc_heads.shape[1]))
    )

    final_elev = (
        np.vstack(all_elev)
        if all_elev
        else np.empty((0, bc_elev.shape[1]))
    )

    return final_parts, final_heads, final_elev


def get_outer_mesh_boundary(msh_outline):
    """
    Extract the exterior boundary of a polygon mesh,
    excluding all interior holes.

    Returns
    -------
    GeoDataFrame
        One closed LineString representing the outer boundary.
        CRS is preserved.
    """

    if msh_outline.crs is None:
        raise ValueError("Input mesh must have a defined CRS.")

    geom = unary_union(msh_outline.geometry.dropna())

    if isinstance(geom, Polygon):
        outer = LineString(geom.exterior.coords)

    elif isinstance(geom, MultiPolygon):
        # Select largest connected polygon
        largest = max(geom.geoms, key=lambda g: g.area)
        outer = LineString(largest.exterior.coords)

    else:
        raise ValueError(
            f"Expected Polygon or MultiPolygon, got {geom.geom_type}"
        )

    return gpd.GeoDataFrame(
        geometry=[outer],
        crs=msh_outline.crs,
    )