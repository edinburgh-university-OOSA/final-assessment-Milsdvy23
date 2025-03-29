#Task_One
'''Include an introduction to explain what the report will do. Use the starting code
to develop a script capable of reading in an LVIS L1B file and plotting a graph of
an arbitrary waveform (any index, provided by you at run time). Include a
description of the code developed and an example command to show how to
produce the graph in your repository’s README, along with an example of the
output.'''

from lvisClass import lvisData
import numpy as np
from pyproj import Transformer
import argparse
import matplotlib.pyplot as plt

def getCmdArgs():
    '''Get command line arguments'''
    parser = argparse.ArgumentParser(description="Plot a specific waveform from an LVIS file")
    parser.add_argument("filename", type=str, help="Path to the LVIS file")
    parser.add_argument("--index", type=int, required=True, help="Index of the waveform to plot")
    return parser.parse_args()

class plotLVIS(lvisData):
    '''A class inheriting from lvisData to add a plotting method'''

    def reprojectLVIS(self, outEPSG):
        '''A method to reproject the geolocation data'''
        reproject = Transformer.from_crs('epsg:4326', outEPSG, always_xy=True)
        self.x, self.y = reproject.transform(self.lon, self.lat)

    def plotWave(self, index):
        '''Method to plot a single waveform'''
        print(f"Plotting waveform at index: {index}")
        plt.figure()
        plt.plot(self.waves[index], self.z[index])
        plt.title(f'Waveform {index}')
        plt.xlabel('Intensity of RXWaveForm')
        plt.ylabel('Elevation (m)')
        plt.grid(True)
        plt.savefig(f'Waveform_{index}.png')
        plt.show()
        print(f"Graph saved as Waveform_{index}.png")

if __name__ == "__main__":
    '''Main block'''
    args = getCmdArgs()
    filename = args.filename
    index = args.index

    b = plotLVIS(filename, onlyBounds=True)

    x0 = b.bounds[0]
    y0 = b.bounds[1]
    x1 = (b.bounds[2] - b.bounds[0]) / 20 + b.bounds[0]
    y1 = (b.bounds[3] - b.bounds[1]) / 20 + b.bounds[1]

    lvis = plotLVIS(filename, minX=x0, minY=y0, maxX=x1, maxY=y1)

    lvis.setElevations()
    lvis.reprojectLVIS('epsg:3031')

    if index < 0 or index >= len(lvis.waves):
        print("Index out of range.")
    else:
        lvis.plotWave(index)

# command to call matplotlib graph: python \#Task_One.py --index 2000