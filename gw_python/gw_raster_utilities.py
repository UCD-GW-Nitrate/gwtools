import numpy as np
from shapely.geometry import mapping
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.mask import mask


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