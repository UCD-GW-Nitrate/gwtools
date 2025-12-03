import warnings
import geopandas as gpd
import numpy as np

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