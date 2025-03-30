#Task_One
from lvisClass import lvisData # Imports file for handling lvis data. 
import numpy as np # used for array manipulation. 
from pyproj import Transformer # used for reprojecting coordinates. 
import argparse # enables command lines arguments.
import matplotlib.pyplot as plt # plots the graphs.

def getCmdArgs():
    # Specifying the command line arguments so nothing is hard coded.
    '''Get command line arguments'''
    parser = argparse.ArgumentParser(description="Plot a specific waveform from an LVIS file")
    parser.add_argument("filename", type=str, help="Path to the LVIS file")
    parser.add_argument("--index", type=int, required=True, help="Index of the waveform to plot")
    return parser.parse_args()

# Utilising object oriented programming
class plotLVIS(lvisData):
    '''A class inheriting from lvisData to add a plotting method usng matplotlib'''

    def reprojectLVIS(self, outEPSG):
        '''A method to reproject the geolocation data'''
        # this converts coordinates from epsg: 4326 (which uses google earth projections) to a specified epsg, in this case 3031 which covers antarctica.
        reproject = Transformer.from_crs('epsg:4326', outEPSG, always_xy=True)
        # This applies to transformation to the x and y coordinates.
        self.x, self.y = reproject.transform(self.lon, self.lat)

    def plotWave(self, index):
        '''Method to plot a single waveform depending on the specified index'''
        # This plots the wave's intensity and elevation using matplotlib at the defined index which the user can choose.
        print(f"Plotting waveform at index: {index}")
        plt.figure()
        # Plots the wave intensity and the elevation.
        plt.plot(self.waves[index], self.z[index])
        # Titling and labelling the graph.
        plt.title(f'Waveform {index}')
        plt.xlabel('Intensity of RXWaveForm')
        plt.ylabel('Elevation (m)')
        # Enables a grid for better readability. 
        plt.grid(True)
        # Saves the figure as a png,called the waveform index.
        plt.savefig(f'Waveform_{index}.png')
        # Show the figure.
        plt.show()
        print(f"Graph saved as Waveform_{index}.png")

# ensures the code runs only if it executes correctly. 
if __name__ == "__main__":
    '''Main block'''
    #Parses command line arguments.
    args = getCmdArgs()
    # Retrieves filename and index from the arguments.
    filename = args.filename
    index = args.index

    # Initialises plotLVIS to retrieve the bounding box of the data. 
    b = plotLVIS(filename, onlyBounds=True)

    # Creates a new bounding box, 1/20th of the original size.
    x0 = b.bounds[0]
    y0 = b.bounds[1]
    x1 = (b.bounds[2] - b.bounds[0]) / 20 + b.bounds[0]
    y1 = (b.bounds[3] - b.bounds[1]) / 20 + b.bounds[1]

    # processes the data within this new bounding box. 
    lvis = plotLVIS(filename, minX=x0, minY=y0, maxX=x1, maxY=y1)

    # Sets the elevation and reprojects. 
    lvis.setElevations()
    lvis.reprojectLVIS('epsg:3031')

    # Checks if the index is within the waveform. 
    if index < 0 or index >= len(lvis.waves):
        print("Index out of range.")
    else:
        # plots the waveform.
        lvis.plotWave(index)

# python 'Task_One.py' '/geos/netdata/oosa/assignment/lvis/2009/ILVIS1B_AQ2009_1020_R1408_052195.h5' --index 2000