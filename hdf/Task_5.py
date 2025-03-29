'''Task_5 
You now have two geotiffs of the elevation 6 years apart. Write a new class or
function to read both DEMs in to RAM (see the handleTiff.py scripts provided for
an example of one way to read geotiffs). Analyse the region where both overlap
and determine the total volume change between the two dates. Report this
change in ice volume.
Add a method to produce a new geotiff showing the change in elevation between
the two dates.
Include a description of the code developed and an example command to show
how to produce your graphs and change map in your repository’s README,
along with a map of the difference, produced in python '''

from pyproj import Transformer
import argparse
import matplotlib.pyplot as plt
import os
import numpy as np
from lvisClass import lvisData
from scipy.ndimage import gaussian_filter1d 
from pyproj import Proj, transform # package for reprojecting data
from osgeo import gdal             # pacage for handling geotiff data
from osgeo import osr              # pacage for handling projection information
import numpy as np
from lvisExample import plotLVIS
from rasterio.transform import from_origin
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.enums import Resampling
import cartopy.crs as ccrs
from rasterio.merge import merge
from glob import glob

from osgeo import gdal

from osgeo import gdal
import numpy as np

from osgeo import gdal, gdalconst
import numpy as np

def read_tiff(filename):
    ds = gdal.Open(filename)
    if not ds:
        raise IOError(f"Unable to open file: {filename}")
    data = ds.GetRasterBand(1).ReadAsArray()
    transform = ds.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    return data, xOrigin, yOrigin, pixelWidth, pixelHeight, ds.RasterXSize, ds.RasterYSize

def resample_to_match(ds, match_ds):
    resampled = gdal.GetDriverByName('MEM').Create(
        '', match_ds.RasterXSize, match_ds.RasterYSize, 1, gdalconst.GDT_Float32)
    resampled.SetGeoTransform(match_ds.GetGeoTransform())
    resampled.SetProjection(match_ds.GetProjection())
    gdal.ReprojectImage(ds, resampled, ds.GetProjection(), match_ds.GetProjection(), gdalconst.GRA_Bilinear)
    return resampled

def calculate_overlap(data1, data2, metadata1, metadata2):
    _, xOrigin1, yOrigin1, pixelWidth1, pixelHeight1, _, _ = metadata1
    _, xOrigin2, yOrigin2, pixelWidth2, pixelHeight2, _, _ = metadata2
    x_start = max(xOrigin1, xOrigin2)
    y_start = min(yOrigin1, yOrigin2)
    x_end = min(xOrigin1 + data1.shape[1] * pixelWidth1, xOrigin2 + data2.shape[1] * pixelWidth2)
    y_end = max(yOrigin1 + data1.shape[0] * pixelHeight1, yOrigin2 + data2.shape[0] * pixelHeight2)
    x_offset1 = int((x_start - xOrigin1) / pixelWidth1)
    y_offset1 = int((yOrigin1 - y_start) / abs(pixelHeight1))
    x_offset2 = int((x_start - xOrigin2) / pixelWidth2)
    y_offset2 = int((yOrigin2 - y_start) / abs(pixelHeight2))
    nX_overlap = int((x_end - x_start) / pixelWidth1)
    nY_overlap = int((y_start - y_end) / abs(pixelHeight1))
    overlap_data1 = data1[y_offset1:y_offset1 + nY_overlap, x_offset1:x_offset1 + nX_overlap]
    overlap_data2 = data2[y_offset2:y_offset2 + nY_overlap, x_offset2:x_offset2 + nX_overlap]
    return overlap_data1, overlap_data2, pixelWidth1, pixelHeight1

def calculate_volume_change(data1, data2, pixelWidth, pixelHeight):
    difference = data2 - data1
    area_per_pixel = abs(pixelWidth * pixelHeight)
    volume_change = np.sum(difference) * area_per_pixel
    return volume_change

def write_tiff(output_filename, data, input_ds):
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_filename, input_ds.RasterXSize, input_ds.RasterYSize, 1, gdalconst.GDT_Float32)
    out_ds.SetGeoTransform(input_ds.GetGeoTransform())
    out_ds.SetProjection(input_ds.GetProjection())
    out_ds.GetRasterBand(1).WriteArray(data)
    out_ds.FlushCache()

import matplotlib.pyplot as plt
import rasterio

filename1 = '/home/s2478921/filled_output2009.tif'
ds1 = gdal.Open(filename1)
data1, xOrigin1, yOrigin1, pixelWidth1, pixelHeight1, nX1, nY1 = read_tiff(filename1)
metadata1 = (data1, xOrigin1, yOrigin1, pixelWidth1, pixelHeight1, nX1, nY1)

filename2 = '/home/s2478921/filled_output2015.tif'
ds2 = gdal.Open(filename2)
resampled_ds2 = resample_to_match(ds2, ds1)
data2 = resampled_ds2.GetRasterBand(1).ReadAsArray()
metadata2 = (data2, xOrigin1, yOrigin1, pixelWidth1, pixelHeight1, nX1, nY1)

overlap_data1, overlap_data2, pixelWidth, pixelHeight = calculate_overlap(data1, data2, metadata1, metadata2)

if overlap_data1.shape == overlap_data2.shape:
    volume_change = calculate_volume_change(overlap_data1, overlap_data2, pixelWidth, pixelHeight)
    print(f"Total Ice Volume Change: {volume_change} cubic meters")
else:
    print("Error: Unable to match the overlap region shapes.")

elevation_difference = data2 - data1

output_filename = '/home/s2478921/elevation_change.tif'
write_tiff(output_filename, elevation_difference, ds1)
print(f"Elevation change GeoTIFF saved as: {output_filename}")
