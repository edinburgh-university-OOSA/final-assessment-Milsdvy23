# import the HDF5 data handler class
 
 from lvisClass import lvisData
 
 ##########################################
 import numpy as np
 from pyproj import Transformer
 import argparse
 import matplotlib.pyplot as plt
 import os
 
 class plotLVIS(lvisData):
   '''A class, ineriting from lvisData
      and add a plotting method'''
 
   def reprojectLVIS(self,outEPSG):
     '''A method to reproject the geolocation data'''
      # set projections
 
     # define projections (epsg needs to be a string)
     reproject = Transformer.from_crs('epsg:4326',outEPSG,always_xy=True)
 
     # reproject data
     x,y=reproject.transform(self.lon,self.lat)
 
     return(x,y)
 
   def plotWaves(self):
 
     '''A method to plot all waveforms'''
 
     for i in range(10):
        plt.figure()
        plt.plot(self.waves[i], self.z[i])
        plt.title(f'Waveform {i+1}')
        plt.xlabel('Intensity of RXWaveForm')
        plt.ylabel('Elevation (m)')
        plt.grid(True)
        plt.show()
   
 
 ##########################################
 
 if __name__=="__main__":
   '''Main block'''
 
   filename='/geos/netdata/oosa/assignment/lvis/2009/ILVIS1B_AQ2009_1020_R1408_049700.h5'
 
   # create instance of class with "onlyBounds" flag
   b=plotLVIS(filename,onlyBounds=True)
 
   # to make a MWE,
   # from the total file bounds
   # choose a spatial subset
   x0=b.bounds[0]
   y0=b.bounds[1]
   x1=(b.bounds[2]-b.bounds[0])/20+b.bounds[0]
   y1=(b.bounds[3]-b.bounds[1])/20+b.bounds[1]
 
 
   # read in all data within our spatial subset
   lvis=plotLVIS(filename,minX=x0,minY=y0,maxX=x1,maxY=y1)
 
   # set elevation
   lvis.setElevations()
 
   # reproject the data
 
   lvis.reprojectLVIS('epsg:3031')
 
   lvis.plotWaves()
 
   #print(lvis.waves)
 
   #print(np.shape(lvis.waves))
 
   # plot up some waveforms using your new method
