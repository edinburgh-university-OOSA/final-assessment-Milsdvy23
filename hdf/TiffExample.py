'''TiffExample'''
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
import rasterio as rasterio
from rasterio.transform import from_origin
from glob import glob

def writeTiff(self,data,filename="chm.tif",epsg=27700):
    '''
    Write a geotiff from a raster layer
    '''
    # set geolocation information (note geotiffs count down from top edge in Y)
    geotransform = (self.minX, self.res, 0, self.maxY, 0, -1*self.res)

    # load data in to geotiff object
    dst_ds = gdal.GetDriverByName('GTiff').Create(filename, self.nX, self.nY, 1, gdal.GDT_Float32)

    dst_ds.SetGeoTransform(geotransform)    # specify coords
    srs = osr.SpatialReference()            # establish encoding
    srs.ImportFromEPSG(epsg)                # WGS84 lat/long
    dst_ds.SetProjection(srs.ExportToWkt()) # export coords to file
    dst_ds.GetRasterBand(1).WriteArray(data)  # write image to the raster
    dst_ds.GetRasterBand(1).SetNoDataValue(-999)  # set no data value
    dst_ds.FlushCache()                     # write to disk
    dst_ds = None

    print("Image written to",filename)
    return 