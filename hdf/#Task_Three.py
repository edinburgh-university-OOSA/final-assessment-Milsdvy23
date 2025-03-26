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
  
    #Make a geotiff from an array of points with average elevation calculation
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
    imageArr[imageArr == -999.0] = np.nan
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
  
    fileList = glob('/geos/netdata/oosa/assignment/lvis/2015/*1.h5')
    output_dir = '/home/s2478921/oose/final-assessment-Milsdvy23/hdf'

    for filename in fileList:

         #'''b=lvisGround(filename,onlyBounds=True)

         #x0=b.bounds[0]
         #y0=b.bounds[1]
         #x1=(b.bounds[2]-b.bounds[0])/1+b.bounds[0]
         #y1=(b.bounds[3]-b.bounds[1])/1+b.bounds[1]

         #lvis=lvisGround(filename,minX=x0,minY=y0,maxX=x1,maxY=y1)

         #lvis.setElevations()

         #lvis.estimateGround()

         #lvis.writeTiff(data=lvis.zG, x=lvis.lon, y=lvis.lat, res=40)

         #base_name = os.path.basename(filename).replace('.h5', '.tif')
         #output_filename = os.path.join(output_dir, base_name)

         #lvis.writeTiff(data=lvis.zG, x=lvis.lon, y=lvis.lat, res=0.001, epsg=4326, filename=output_filename)'''

         tif = glob('/home/s2478921/oose/final-assessment-Milsdvy23/hdf/*1.tif')

import os
from glob import glob
from osgeo import gdal

# Define paths
dir_path = '/home/s2478921/oose/final-assessment-Milsdvy23/hdf/'
out_path = '/home/s2478921/scratch/mosaic3.tif'

# Get list of TIFF files
file_list = glob(os.path.join(dir_path, '*1.tif'))

# Check if files were found
if not file_list:
    print("No TIFF files found.")
else:
    print(f"Found {len(file_list)} TIFF files.")

# Use GDAL to merge the files
vrt = gdal.BuildVRT('merged.vrt', file_list)

if vrt is not None:
    gdal.Translate(out_path, vrt, format='GTiff')
    print(f"Mosaic saved to {out_path}")
else:
    print("Failed to create VRT.")

# Clean up
vrt = None



'''import os
import numpy as np
import glob
import rasterio
from rasterio.merge import merge
from rasterio.plot import show

# Define input and output paths
dir_path = '/home/s2478921/oose/final-assessment-Milsdvy23/hdf/'
out_path = '/home/s2478921/scratch/mosaic2.tif'

# Define search criteria for TIF files
search_criteria = "*1.tif"
q = os.path.join(dir_path, search_criteria)
print(q)
# Get list of file paths matching the search criteria
dem_fps = glob.glob(q)

print(dem_fps)

src_files_to_mosaic = []

for fp in dem_fps:
  src = rasterio.open(fp)
  src_files_to_mosaic.append(src)

print(src_files_to_mosaic)

mosaic, out_trans = merge(src_files_to_mosaic)

plt.imshow(mosaic[0], cmap='terrain')
plt.show()

out_meta = src.meta.copy()

out_meta.update({"driver": "GTiff",
                     "height": mosaic.shape[1],
                   "width": mosaic.shape[2],
                     "transform": out_trans,
                     "count": mosaic.shape[0],
                     "dtype": mosaic.dtype,
                     }
                    )

with rasterio.open(out_path, "w", **out_meta) as dest:
  dest.write(mosaic)'''










