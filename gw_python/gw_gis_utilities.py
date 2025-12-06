import pandas as pd
import geopandas as gpd

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