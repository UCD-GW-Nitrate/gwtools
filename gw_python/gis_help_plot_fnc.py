import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon, LineString, MultiLineString
from shapely.geometry import box

def plot_shape_and_nodes(gdf, i, ax=None):
    """
    Plot the geometry of gdf at index i together with its vertices.
    Labels each vertex with its sequential order.
    """

    geom = gdf.loc[i, "geometry"]

    # Create axis if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    # --- Plot the geometry ----------------------------------------------------
    gpd.GeoSeries([geom]).plot(ax=ax, color='lightgray', edgecolor='black')

    # --- Extract coordinates depending on geometry type -----------------------
    coords = []

    if isinstance(geom, Polygon):
        coords = list(geom.exterior.coords)

    elif isinstance(geom, MultiPolygon):
        # take the largest polygon
        largest = max(geom.geoms, key=lambda p: p.area)
        coords = list(largest.exterior.coords)

    elif isinstance(geom, LineString):
        coords = list(geom.coords)

    elif isinstance(geom, MultiLineString):
        # take the longest line
        longest = max(geom.geoms, key=lambda l: l.length)
        coords = list(longest.coords)

    else:
        raise TypeError(f"Unsupported geometry type: {type(geom)}")

    # --- Plot points and indices ---------------------------------------------
    xs, ys = zip(*coords)
    ax.scatter(xs, ys, s=20, color='red')

    for idx, (x, y) in enumerate(coords):
        ax.text(x, y, str(idx), fontsize=12, color='blue')

    ax.set_title(f"Geometry and Vertex Order for Row {i}")
    ax.set_aspect('equal')
    return ax

def plot_streamline_and_overlapping_mesh(S_gdf, mesh_gdf, idx, ax=None):
    selected_lines = S_gdf.iloc[idx]

    minx, miny, maxx, maxy = selected_lines.total_bounds
    bbox = box(minx, miny, maxx, maxy)

    overlapping_polygons = mesh_gdf[mesh_gdf.intersects(bbox)]

    overlapping_polygons.plot(ax=ax, facecolor="none", edgecolor="orange", linewidth=0.5)

    selected_lines.plot(ax=ax, color="blue", linewidth=2)

    #gpd.GeoSeries([bbox]).plot(ax=ax, facecolor="none", edgecolor="red", linestyle="--")

    #ax.set_title("Selected LineStrings and Overlapping Mesh Polygons", fontsize=14)
    ax.set_aspect("equal")