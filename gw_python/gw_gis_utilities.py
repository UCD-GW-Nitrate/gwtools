import pandas as pd
import geopandas as gpd
import numpy as np
from scipy.spatial import Delaunay
from pathlib import Path

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