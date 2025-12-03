import warnings
import geopandas as gpd

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