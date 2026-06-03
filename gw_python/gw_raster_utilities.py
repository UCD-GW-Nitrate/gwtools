import numpy as np
from shapely.geometry import mapping
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.mask import mask
from rasterio.features import rasterize
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
    template_raster,
    out_raster,
    nodata=-9999.0,
    all_touched=False,
    dtype="float32",
):
    """
    Rasterize one GeoDataFrame field using the grid, transform, and CRS
    of a template raster.

    The GeoDataFrame and template raster are expected to use the same CRS.
    """

    template_raster = Path(template_raster)
    out_raster = Path(out_raster)

    with rasterio.open(template_raster) as src:
        meta = src.meta.copy()
        transform = src.transform
        out_shape = (src.height, src.width)
        template_crs = src.crs

    if gdf.crs != template_crs:
        raise ValueError(
            f"CRS mismatch: gdf.crs={gdf.crs}, template_crs={template_crs}"
        )

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

    meta.update({
        "driver": "GTiff",
        "count": 1,
        "dtype": dtype,
        "nodata": nodata,
        "compress": "lzw",
    })

    out_raster.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.open(out_raster, "w", **meta) as dst:
        dst.write(arr, 1)

    return arr