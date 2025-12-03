import geopandas as gpd
from shapely.geometry import Polygon

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

