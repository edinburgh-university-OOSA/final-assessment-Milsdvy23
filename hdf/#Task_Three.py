#Task_Three
'''
Adapt the code from task 1 to process all of the 2015 data in to a single DEM, in
geotiff format, at a resolution of your choice. Choose an appropriate resolution to
keep RAM usage and the output file size below 2 Gbytes. Describe how you chose
to combine the files and why you chose the final resolution. This may require
subsetting to keep the RAM usage below 2 Gbytes (do not go over 2 Gbytes)
and may need to be run overnight. Your code will be assessed and run to see
whether you stayed beneath this limit. Include a description of the tests you
performed on your code before running over the full dataset. If processing is
likely to take too long (longer than overnight), it is acceptable to choose a spatial
subset small enough that processing can be completed overnight, and only
process that for your report. If you do subset, make sure that your subset will
include data from both 2009 and 2015 so that you can answer task 4.
Include a description of the code developed and an example command to show
how to produce your DEM in your repository’s README, along with images of
your DEM '''

#######################################

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

#######################################

class lvisGround(lvisData):
  '''
  LVIS class with extra processing steps
  to allow it to found the ground over ice
  '''
  #######################################################

  def estimateGround(self,sigThresh=5,statsLen=10,minWidth=3,sWidth=0.5):
    '''
    Processes waveforms to estimate ground
    Only works for bare Earth. DO NOT USE IN TREES
    '''

    # find noise statistics
    self.findStats(statsLen=statsLen)

    # set threshold
    threshold=self.setThreshold(sigThresh)

    # remove background
    self.denoise(threshold,minWidth=minWidth,sWidth=sWidth)

    # find centre of gravity of remaining signal
    self.CofG()

  #######################################################

  def setThreshold(self,sigThresh):
    '''
    Set a noise threshold
    '''
    threshold=self.meanNoise+sigThresh*self.stdevNoise
    return(threshold)

  #######################################################

  def CofG(self):
    '''
    Find centre of gravity of denoised waveforms
    sets this to an array of ground elevation
    estimates, zG
    '''

    # allocate space and put no data flags
    self.zG=np.full((self.nWaves),-999.0)

    # loop over waveforms
    for i in range(0,self.nWaves):
      if(np.sum(self.denoised[i])>0.0):   # avoid empty waveforms (clouds etc)
        self.zG[i]=np.average(self.z[i],weights=self.denoised[i])  # centre of gravity

  ##############################################

  def findStats(self,statsLen=10):
    '''
    Finds standard deviation and mean of noise
    '''
    # make empty arrays
    self.meanNoise=np.empty(self.nWaves)
    self.stdevNoise=np.empty(self.nWaves)

    # determine number of bins to calculate stats over
    res=(self.z[0,0]-self.z[0,-1])/self.nBins    # range resolution
    noiseBins=int(statsLen/res)   # number of bins within "statsLen"

    # loop over waveforms
    for i in range(0,self.nWaves):
      self.meanNoise[i]=np.mean(self.waves[i,0:noiseBins])
      self.stdevNoise[i]=np.std(self.waves[i,0:noiseBins])

  ##############################################

  def denoise(self,threshold,sWidth=0.5,minWidth=3):
    '''
    Denoise waveform data
    '''
    res=(self.z[0,0]-self.z[0,-1])/self.nBins  

    self.denoised=np.full((self.nWaves,self.nBins),0)

    for i in range(0,self.nWaves):
      print("Denoising wave",i+1,"of",self.nWaves)

      self.denoised[i]=self.waves[i]-self.meanNoise[i]

      self.denoised[i,self.denoised[i]<threshold[i]]=0.0

      # minimum acceptable width
      binList=np.where(self.denoised[i]>0.0)[0]
      for j in range(0,binList.shape[0]):       # loop over waveforms
        if((j>0)&(j<(binList.shape[0]-1))):    # are we in the middle of the array?
          if((binList[j]!=binList[j-1]+1)|(binList[j]!=binList[j+1]-1)):  # are the bins consecutive?
            self.denoised[i,binList[j]]=0.0   # if not, set to zero

      # smooth
      self.denoised[i]=gaussian_filter1d(self.denoised[i],sWidth/res)

  def writeTiff(self, data, x, y, res, filename="lvis_image.tif", epsg=4326):

    minX = np.min(x)
    maxX = np.max(x)
    minY = np.min(y)
    maxY = np.max(y)

    nX = int((maxX - minX) / res + 1)
    nY = int((maxY - minY) / res + 1)

    imageArr = np.full((nY, nX), -999.0)
    elevationSum = np.full((nY, nX), 0.0)
    count = np.full((nY, nX), 0)

    xInds = np.array(np.floor((x - minX) / res), dtype=int)
    yInds = np.array(np.floor((maxY - y) / res), dtype=int)

    for i in range(len(data)):
        elevationSum[yInds[i], xInds[i]] += data[i]
        count[yInds[i], xInds[i]] += 1

    validPixels = count > 0
    imageArr[validPixels] = elevationSum[validPixels] / count[validPixels]
    imageArr[imageArr == -999.0]
    geotransform = (minX, res, 0, maxY, 0, -res)

    dst_ds = gdal.GetDriverByName('GTiff').Create(filename, nX, nY, 1, gdal.GDT_Float32)
    dst_ds.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    dst_ds.SetProjection(srs.ExportToWkt())
    dst_ds.GetRasterBand(1).WriteArray(imageArr)
    dst_ds.GetRasterBand(1).SetNoDataValue(-999)
    dst_ds.FlushCache()
    dst_ds = None

    print("Image written to", filename)
    return
  
def to_360_range(lon):
    return lon if lon >= 0 else lon + 360

import gc

def process_files_in_batches(file_list, batch_size, scratch_dir):

    for i in range(0, len(file_list), batch_size):
        batch = file_list[i:i + batch_size]
        for filename in batch:
            b = lvisGround(filename, onlyBounds=True)

            file_min_lon = to_360_range(b.bounds[0])
            file_max_lon = to_360_range(b.bounds[2])
            file_min_lat = b.bounds[1]
            file_max_lat = b.bounds[3]

            overlaps_lon = (file_min_lon < fixed_max_lon and file_max_lon > fixed_min_lon)
            overlaps_lat = (file_min_lat < fixed_max_lat and file_max_lat > fixed_min_lat)

            if overlaps_lon and overlaps_lat:
                lvis = lvisGround(
                   filename, minX=file_min_lon, minY=file_min_lat,
                   maxX=file_max_lon, maxY=file_max_lat)

                lvis.setElevations()
                lvis.estimateGround()

                lvis.writeTiff(data=lvis.zG, x=lvis.lon, y=lvis.lat, res=150)

                base_name = os.path.basename(filename).replace('.h5', '.tif')
                output_filename = os.path.join(scratch_dir, base_name)

                lvis.writeTiff(data=lvis.zG, x=lvis.lon, y=lvis.lat, res=150, epsg=4326, filename=output_filename)

                del lvis
                gc.collect()

import rasterio
from rasterio.merge import merge
from glob import glob
import os

def merge_tiff_files(scratch_dir, output_tiff_path):
    
    file_list = glob(os.path.join(scratch_dir, "*.tif"))
    
    src_files_to_mosaic = [rasterio.open(fp) for fp in file_list]

  
    mosaic, out_trans = merge(src_files_to_mosaic)
    
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "crs": src_files_to_mosaic[0].crs
    })

    with rasterio.open(output_tiff_path, "w", **out_meta) as dest:
        dest.write(mosaic)
    
    for src in src_files_to_mosaic:
        src.close()

from scipy.ndimage import convolve

def fill_gaps_in_dem(input_path, output_path, nodata_value=-999):
    
    with rasterio.open(input_path) as src:
        data = src.read(1)
        mask = data == nodata_value

        kernel = np.ones((3, 3))

        mask_inv = np.logical_not(mask)
        convolved = convolve(data * mask_inv, kernel, mode='constant', cval=0.0)
        weight = convolve(mask_inv.astype(float), kernel, mode='constant', cval=0.0)
        
        filled_data = np.where(mask, convolved / np.maximum(weight, 1), data)

        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "dtype": 'float32',
            "nodata": nodata_value
        })

        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(filled_data, 1)

import psutil

def check_memory():
    process = psutil.Process()
    mem_info = process.memory_info()
    return mem_info.rss / (1024 ** 2)

if __name__ == "__main__":
    fileList = glob('/geos/netdata/oosa/assignment/lvis/2015/*.h5')
    scratch_dir = '/home/s2478921/scratch/'
    out_path = '/home/s2478921/scratch/merged_output.tif'
    filled_output_path = '/home/s2478921/filled_output.tif'
    batch_size = 2

    center_lat = -74.8976
    center_lon = -99.9302
    east_margin = 0.5
    west_margin = 0.5
    north_south_margin = 0.5
    fixed_min_lon = to_360_range(center_lon - west_margin)
    fixed_max_lon = to_360_range(center_lon + east_margin)
    fixed_min_lat = center_lat - north_south_margin
    fixed_max_lat = center_lat + north_south_margin

    print(f"Initial Memory Usage: {check_memory()} MB")

    #process_files_in_batches(fileList, batch_size, scratch_dir)

    #merge_tiff_files('/home/s2478921/scratch', '/home/s2478921/merged.vrt')
    #fill_gaps_in_dem('/home/s2478921/merged.vrt', '/home/s2478921/filled_output.tif')

    #plot_dem(filled_output_path)



  











