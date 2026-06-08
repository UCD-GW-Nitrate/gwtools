import numpy as np
from shapely.geometry import Polygon
from shapely import affinity
import geopandas as gpd

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


def clip_elements_to_rectangle(elem_gdf, rect_poly):
    """
    Clip element polygons to a rectangular model domain and identify
    which elements were cut by the boundary.

    Parameters
    ----------
    elem_gdf : GeoDataFrame
        Input element polygons.
    rect_poly : shapely.geometry.Polygon
        Rectangle polygon defining the local model domain.

    Returns
    -------
    clipped_gdf : GeoDataFrame
        Clipped elements with the following added fields:

        - is_internal : bool
            Original element was completely inside the rectangle.
        - is_boundary_clipped : bool
            Original element intersected the rectangle boundary and
            was clipped.
    """

    # Identify element status before clipping
    inside = elem_gdf.within(rect_poly)
    intersect = elem_gdf.intersects(rect_poly)

    # Keep only intersecting elements
    work_gdf = elem_gdf.loc[intersect].copy()

    # Add flags
    work_gdf["is_internal"] = inside.loc[intersect].values
    work_gdf["is_boundary_clipped"] = (~inside.loc[intersect]).values

    # Clip geometries
    rect_gdf = gpd.GeoDataFrame(
        geometry=[rect_poly],
        crs=elem_gdf.crs
    )

    clipped_gdf = gpd.clip(work_gdf, rect_gdf)

    return clipped_gdf.reset_index(drop=True)