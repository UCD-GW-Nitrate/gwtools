import numpy as np
from shapely.geometry import mapping
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.mask import mask
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.crs import CRS
from pathlib import Path


def add_raster_values_to_df(df, raster_path, x_col='X', y_col='Y', field_name='initConc'):
    """
    Interpolates raster values at given coordinates and appends them to the DataFrame.

    Parameters:
        df (pd.DataFrame): DataFrame containing X and Y coordinate columns.
        raster_path (str): Path to the raster file.
        x_col (str): Name of the column with X coordinates.
        y_col (str): Name of the column with Y coordinates.
        field_name (str): Name of the new column to store raster values.

    Returns:
        pd.DataFrame: Updated DataFrame with the new column.
    """
    coords = list(zip(df[x_col], df[y_col]))

    with rasterio.open(raster_path) as src:
        # Use WarpedVRT for potential resampling/interpolation
        with WarpedVRT(src, resampling=Resampling.bilinear) as vrt:
            values = [
                list(vrt.sample([(x, y)]))[0][0]  # sample returns a generator of 1D arrays
                for x, y in coords
            ]

    df[field_name] = values
    return df

def clip_raster_to_poly(raster_file_name, poly):
    geoms = [mapping(geom) for geom in poly.geometry]

    with rasterio.open(raster_file_name) as src:
        A_raster_clipped, A_raster_transform = mask(
            src, geoms, crop=True, all_touched=False, filled=True, nodata=np.nan
        )
        R_raster_clipped = src.profile.copy()
        R_raster_clipped.update({
            "height": A_raster_clipped.shape[1],
            "width": A_raster_clipped.shape[2],
            "transform": A_raster_transform
        })

    return A_raster_clipped, R_raster_clipped

def write_raster_like_ref(data, filename, ref_raster):
    with rasterio.open(ref_raster) as ref:
        meta = ref.meta.copy()
        ref_height, ref_width = ref.height, ref.width

    # Validate dimensions
    if data.shape != (ref_height, ref_width):
        raise ValueError(
            f"Data shape {data.shape} does not match reference raster "
            f"shape ({ref_height}, {ref_width})."
        )

    # Update metadata
    meta.update({
        "count": 1,
        "height": ref_height,
        "width": ref_width,
        "dtype": data.dtype,
        "nodata": np.nan
    })

# Write output GeoTIFF
    with rasterio.open(filename, "w", **meta) as dst:
        dst.write(data, 1)

def rasterize_gdf_field(
    gdf,
    field,
    grid_info,
    out_raster,
    nodata=-9999.0,
    all_touched=False,
    dtype="float32",
):
    """
    Rasterize one GeoDataFrame field using grid information from a dictionary.

    Expected grid_info keys:
        x_origin, y_origin, cellsize_x, cellsize_y, nrows, ncols

    The output raster CRS is taken from gdf.crs.
    """

    out_raster = Path(out_raster)

    if gdf.crs is None:
        raise ValueError("Input GeoDataFrame has no CRS defined.")

    x_origin = grid_info["x_origin"]
    y_origin = grid_info["y_origin"]
    cellsize_x = grid_info["cellsize_x"]
    cellsize_y = grid_info["cellsize_y"]
    nrows = grid_info["nrows"]
    ncols = grid_info["ncols"]

    transform = from_origin(
        x_origin,
        y_origin,
        cellsize_x,
        cellsize_y,
    )

    out_shape = (nrows, ncols)

    shapes = (
        (geom, value)
        for geom, value in zip(gdf.geometry, gdf[field])
        if geom is not None and not geom.is_empty and np.isfinite(value)
    )

    arr = rasterize(
        shapes=shapes,
        out_shape=out_shape,
        transform=transform,
        fill=nodata,
        dtype=dtype,
        all_touched=all_touched,
    )

    meta = {
        "driver": "GTiff",
        "height": nrows,
        "width": ncols,
        "count": 1,
        "dtype": dtype,
        "crs": gdf.crs,
        "transform": transform,
        "nodata": nodata,
        "compress": "lzw",
    }

    out_raster.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.open(out_raster, "w", **meta) as dst:
        dst.write(arr, 1)

    return arr

def get_raster_grid_info(raster_file):
    """
    Read basic grid information from a raster.

    Returns
    -------
    dict
        {
            "x_origin": float,
            "y_origin": float,
            "cellsize_x": float,
            "cellsize_y": float,
            "ncols": int,
            "nrows": int,
            "xmin": float,
            "ymin": float,
            "xmax": float,
            "ymax": float,
            "crs": rasterio.crs.CRS
        }
    """
    raster_file = Path(raster_file)

    with rasterio.open(raster_file) as src:
        transform = src.transform
        bounds = src.bounds

        return {
            "x_origin": transform.c,          # upper-left x
            "y_origin": transform.f,          # upper-left y
            "cellsize_x": transform.a,
            "cellsize_y": abs(transform.e),
            "ncols": src.width,
            "nrows": src.height,
            "xmin": bounds.left,
            "ymin": bounds.bottom,
            "xmax": bounds.right,
            "ymax": bounds.top,
            "crs": src.crs,
        }

def read_raster_array(raster_file, masked=False):
    """
    Read a raster and return it as a NumPy array.

    Parameters
    ----------
    raster_file : str or Path
    masked : bool, default=False
        If True, returns a masked array using the raster nodata value.

    Returns
    -------
    np.ndarray or np.ma.MaskedArray
    """
    raster_file = Path(raster_file)

    with rasterio.open(raster_file) as src:
        arr = src.read(1, masked=masked)

    return arr

def write_array_like_raster(
    arr,
    reference_raster,
    out_raster,
    nodata=-9999.0,
    dtype="float32",
    compress="lzw",
    force_crs="EPSG:3310",
):
    reference_raster = Path(reference_raster)
    out_raster = Path(out_raster)

    with rasterio.open(reference_raster) as src:
        transform = src.transform
        height = src.height
        width = src.width
        profile = src.profile.copy()

    if arr.shape != (height, width):
        raise ValueError(
            f"Array shape {arr.shape} does not match reference raster shape "
            f"({height}, {width})"
        )

    crs = CRS.from_user_input(force_crs)

    profile.update({
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": 1,
        "dtype": dtype,
        "crs": crs,
        "transform": transform,
        "nodata": nodata,
        "compress": compress,
    })

    out_raster.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.open(out_raster, "w", **profile) as dst:
        dst.write(arr.astype(dtype), 1)