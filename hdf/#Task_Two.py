#Task_Two
'''Then add the methods needed to process that same LVIS file in to a DEM of any
chosen resolution. This will require a command line parser to provide the
filename and the resolution at runtime, as well as some additional methods
adding to the scripts provided. Make use of objects and inheritance where
possible to make the code easier to maintain. Do not have duplicate code in
your repository.
Use this code to process a flight line of your choice over PIG to a 30 m resolution
DEM in geotiff format. Include the DEM as a figure and a short discussion of any
interesting features in your report. Note that the geotiff example code given uses
data from only a single lidar spot in each geotiff pixel. A more accurate answer
would be to calculate the average elevation of all intersecting lidar shots. This
can be done for extra credit.
Include a description of the code developed and an example command to show
how to produce your DEM in your repository’s README.'''


'''
Some example functions for processing LVIS data
'''

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

    # find resolution
    res=(self.z[0,0]-self.z[0,-1])/self.nBins    # range resolution

    # make array for output
    self.denoised=np.full((self.nWaves,self.nBins),0)

    # loop over waves
    for i in range(0,self.nWaves):
      print("Denoising wave",i+1,"of",self.nWaves)

      # subtract mean background noise
      self.denoised[i]=self.waves[i]-self.meanNoise[i]

      # set all values less than threshold to zero
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
    '''
    Make a geotiff from an array of points with average elevation calculation
    '''
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

if __name__=="__main__":

  filename='/geos/netdata/oosa/assignment/lvis/2009/ILVIS1B_AQ2009_1020_R1408_049700.h5'

  b=lvisGround(filename,onlyBounds=True)

  x0=b.bounds[0]
  y0=b.bounds[1]
  x1=(b.bounds[2]-b.bounds[0])/20+b.bounds[0]
  y1=(b.bounds[3]-b.bounds[1])/20+b.bounds[1]

  lvis=lvisGround(filename,minX=x0,minY=y0,maxX=x1,maxY=y1)

  lvis.setElevations()

  lvis.estimateGround()

  lvis.writeTiff(data=lvis.zG, x=lvis.lon, y=lvis.lat, res=0.01)
  lvis.writeRasterio(data=lvis.zG, x=lvis.lon, y=lvis.lat, res=30)
#############################################################