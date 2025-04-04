#Task_Two

#######################################

# Coordinate and Spatial Manipulation
from pyproj import Transformer, Proj, transform
from osgeo import gdal, osr

# Data Handling and Processing
import numpy as np

# Plotting and Visualization
import matplotlib.pyplot as plt

# Command-Line Arguments
import argparse

# LVIS Data Management
from lvisClass import lvisData
from lvisExample import plotLVIS

# Image Processing
from scipy.ndimage import gaussian_filter1d

# File and Raster Handling
import os
import rasterio
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

 
  def writeTiff(self, data,x,y,res,filename="lvis_image.tif",epsg=4326):
      '''
      Make a geotiff from an array of points
      '''
      # determines bounds.
      minX=np.min(x)
      maxX=np.max(x)
      minY=np.min(y)
      maxY=np.max(y)

      # determines image size.
      nX=int((maxX-minX)/res+1)
      nY=int((maxY-minY)/res+1)

      # pack in to array
      imageArr=np.full((nY,nX),-999.0)        # make an array of missing data flags

      # calculate the raster pixel index in x and y
      xInds=np.array(np.floor((x-np.min(x))/res),dtype=int)   # need to force to int type
      yInds=np.array(np.floor((np.max(y)-y)/res),dtype=int)
      # floor rounds down. y is from top to bottom

      # this is a simple pack which will assign a single footprint to each pixel
      imageArr[yInds,xInds]=data

      # set geolocation information (note geotiffs count down from top edge in Y)
      geotransform = (minX, res, 0, maxY, 0, -res)

      # load data in to geotiff object
      dst_ds = gdal.GetDriverByName('GTiff').Create(filename, nX, nY, 1, gdal.GDT_Float32)

      dst_ds.SetGeoTransform(geotransform)    # specify coords
      srs = osr.SpatialReference()            # establish encoding
      srs.ImportFromEPSG(epsg)                # WGS84 lat/long
      dst_ds.SetProjection(srs.ExportToWkt()) # export coords to file
      dst_ds.GetRasterBand(1).WriteArray(imageArr)  # write image to the raster
      dst_ds.GetRasterBand(1).SetNoDataValue(-999)  # set no data value
      dst_ds.FlushCache()                     # write to disk
      dst_ds = None

      print("Image written to",filename)
      return
  
def plotTiff(filename="lvis_image.tif"):
      ''' Plotting tiff as a matplotlib graph using viridis colouring'''
      
      # Reading the first wave, and applying a mask to any non-valid values
      with rasterio.open(filename) as src:
        data = src.read(1, masked=True)
        # plotting a figure size as 10,10.
        plt.figure(figsize=(10, 10))
        # setting colour as viridis.
        plt.imshow(data, cmap='viridis')
        # setting colourbar (legend) based on the elevation.
        plt.colorbar(label='Elevation')
        # titling and labelling using matplotlib.
        plt.title('Raster Image')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.gca().invert_yaxis()
        plt.show()

# only run if the code runs correctly.
if __name__ == "__main__":

# argument parsing
  parser = argparse.ArgumentParser(description='Process LVIS data into a DEM.')
  parser.add_argument('filename', type=str, help='Path to the LVIS file')
  parser.add_argument('resolution', type=float, help='Resolution for the DEM')
  parser.add_argument('output', type=str, help='Output filename for the DEM')
  args = parser.parse_args()

 #Retrieves filename and index from the arguments.
  filename = args.filename
  resolution = args.resolution
  output = args.output

  # Initalises lvisground to draw from filename and the bounds from lvisground.
  b=lvisGround(filename,onlyBounds=True)

  # Creates a new bounding box, 1/20th of the original size.
  x0=b.bounds[0]
  y0=b.bounds[1]
  x1=(b.bounds[2]-b.bounds[0])/20+b.bounds[0]
  y1=(b.bounds[3]-b.bounds[1])/20+b.bounds[1]

  # processes the data within this new bounding box. 
  lvis=lvisGround(filename,minX=x0,minY=y0,maxX=x1,maxY=y1)
 
  # drawing on set elevations function
  lvis.setElevations()

  # drawing on estimate ground function.
  lvis.estimateGround()
  # writing tiff, drawing on elevation, coordinates and stating resolution.
  lvis.writeTiff(data=lvis.zG, x=lvis.lon, y=lvis.lat, res=0.01)

  plotTiff('lvis_image.tif')

#############################################################
