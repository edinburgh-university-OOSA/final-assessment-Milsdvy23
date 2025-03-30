
#######################################
import os
from glob import glob

# Data processing and manipulation
import numpy as np
from scipy.ndimage import gaussian_filter1d

# Geospatial libraries
from pyproj import Transformer, Proj, transform
from osgeo import gdal, osr

# Raster handling
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import calculate_default_transform, reproject
from rasterio.enums import Resampling
from rasterio.merge import merge

# Plotting
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Custom classes
from lvisClass import lvisData
from lvisExample import plotLVIS

# Argument parsing
import argparse

import numpy as np
import rasterio
from rasterio.fill import fillnodata
import matplotlib.pyplot as plt

#######################################

# Importing object orientated software engineering.
class lvisGround(lvisData):
    def estimateGround(self, sigThresh=5, statsLen=10, minWidth=3, sWidth=0.5):
        self.findStats(statsLen=statsLen)
        threshold = self.setThreshold(sigThresh)
        self.denoise(threshold, minWidth=minWidth, sWidth=sWidth)
        self.CofG()

    def setThreshold(self, sigThresh):
        threshold = self.meanNoise + sigThresh * self.stdevNoise
        return threshold

    def CofG(self):
        self.zG = np.full((self.nWaves), -999.0)
        for i in range(0, self.nWaves):
            if np.sum(self.denoised[i]) > 0.0:
                self.zG[i] = np.average(self.z[i], weights=self.denoised[i])

    def findStats(self, statsLen=10):
        self.meanNoise = np.empty(self.nWaves)
        self.stdevNoise = np.empty(self.nWaves)
        res = (self.z[0, 0] - self.z[0, -1]) / self.nBins
        noiseBins = int(statsLen / res)
        for i in range(0, self.nWaves):
            self.meanNoise[i] = np.mean(self.waves[i, 0:noiseBins])
            self.stdevNoise[i] = np.std(self.waves[i, 0:noiseBins])

    def denoise(self, threshold, sWidth=0.5, minWidth=3):
        res = (self.z[0, 0] - self.z[0, -1]) / self.nBins
        self.denoised = np.full((self.nWaves, self.nBins), 0)
        for i in range(0, self.nWaves):
            print("Denoising wave", i + 1, "of", self.nWaves)
            self.denoised[i] = self.waves[i] - self.meanNoise[i]
            self.denoised[i, self.denoised[i] < threshold[i]] = 0.0
            binList = np.where(self.denoised[i] > 0.0)[0]
            for j in range(0, binList.shape[0]):
                if j > 0 and j < (binList.shape[0] - 1):
                    if (binList[j] != binList[j - 1] + 1 or
                            binList[j] != binList[j + 1] - 1):
                        self.denoised[i, binList[j]] = 0.0
            self.denoised[i] = gaussian_filter1d(self.denoised[i], sWidth / res)

    def writeTiff(self, data, x, y, res, filename, epsg=4326):
        minX = np.min(x)
        maxX = np.max(x)
        minY = np.min(y)
        maxY = np.max(y)
        nX = int((maxX - minX) / res + 1)
        nY = int((maxY - minY) / res + 1)
        imageArr = np.full((nY, nX), -999.0)
        xInds = np.array(np.floor((x - np.min(x)) / res), dtype=int)
        yInds = np.array(np.floor((np.max(y) - y) / res), dtype=int)
        imageArr[yInds, xInds] = data
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

# converting longitude to 360 to ensure that graph plots correctly.
def convert_longitude(lon):
    return lon + 360 if lon < 0 else lon

# Processing files (taking input directory and output directory as arguments.)
def process_files(input_dir, output_dir):
    ''' This processes the h5 files of either 2015 or 2009. It calcuates bounds, subdivides into tiles and writes tiffs for the valid files.'''
    # looking for H5 files in input directory.
    fileList = glob(os.path.join(input_dir, '*.h5'))
    # if output directory does not exists create output directory.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    #creating a margin of two
    margin = 2
    # coordinates of pine glacier (converting longitude)
    mid_x = convert_longitude(-100.000)
    mid_y = -75.167
    # shifting the bounds right and down
    shift_right = 2
    shift_down = 0.2
    # calcuating min and max x and y by calcuating mid point plus margin. 
    minX = mid_x - margin + shift_right
    maxX = mid_x + margin + shift_right
    minY = mid_y - margin - shift_down
    maxY = mid_y + margin - shift_down
    # setting tile size. 
    tileSize = 2
    xTiles = int((maxX - minX) / tileSize)
    yTiles = int((maxY - minY) / tileSize)
    # for every 3rd file in filelist (to save memory)
    for index, filePath in enumerate(fileList):
        if index % 3 != 0:
            continue

        # Initalises lvisground to draw from filename and the bounds from lvisground.
        b = lvisGround(filePath, onlyBounds=True)

        # if it is within the bounds.
        if b.bounds[0] > maxX or b.bounds[2] < minX or b.bounds[1] > maxY or b.bounds[3] < minY:
            continue
        # Running through each tile.
        for i in range(xTiles):
            for j in range(yTiles):
                tileMinX = minX + i * tileSize
                tileMaxX = tileMinX + tileSize
                tileMinY = minY + j * tileSize
                tileMaxY = tileMinY + tileSize

                # if the tile is within the bounds.
                if b.bounds[0] > tileMaxX or b.bounds[2] < tileMinX or b.bounds[1] > tileMaxY or b.bounds[3] < tileMinY:
                    continue
                
                # calcuating the tile bounds, drawing from lvisGround.
                lvis = lvisGround(filePath, minX=tileMinX, minY=tileMinY, maxX=tileMaxX, maxY=tileMaxY)

                # only processes if there are available waveforms to process.
                if lvis.nWaves == 0:
                    continue
                
                # setting elevations
                lvis.setElevations()
                #estimating ground.
                lvis.estimateGround()
                # using the original file name to use in output file naming.  
                base_name = os.path.splitext(os.path.basename(filePath))[0]
                # constructing the new filename. 
                output_filename = os.path.join(output_dir, f'{base_name}_tile_{i}_{j}.tif')
                # writing tiff based on elevation. x and y coordinates and a high resolution.
                lvis.writeTiff(data=lvis.zG, x=lvis.lon, y=lvis.lat, res=0.01, filename=output_filename)

def merge_tif(scratch_dir, out_path, year):
    ''' This merges the tiff files, methodology based on this website: https://automating-gis-processes.github.io/CSC18/lessons/L6/raster-mosaic.html '''
    # Looking for tif files in directory. 
    dem_fps = glob(os.path.join(scratch_dir, '*.tif'))
    # combining multiple tiff files into once mosaic using rasterio.
    src_files_to_mosaic = [rasterio.open(fp) for fp in dem_fps]
    mosaic, out_trans = merge(src_files_to_mosaic)
    for src in src_files_to_mosaic:
        src.close()
    # plotting and saving figure using matplotlib. 
    plt.imshow(mosaic[0], cmap='terrain')
    # saving figure named depending on the year. 
    plt.savefig(f'Merged_Tiffs_{year}.png')
    plt.close()
    # saving meta data to ensure the tiffs are plotted correctly. 
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "count": mosaic.shape[0],
        "dtype": mosaic.dtype,
    })
    # saving mosaic to out path with meta data. 
    with rasterio.open(out_path, "w", **out_meta) as dest:
        dest.write(mosaic)

def reproject_to_3031(src_path, dest_path, year):
    ''' This function reprojects to 3031 (Antarctica)'''
    # opening the source raster file. 
    with rasterio.open(src_path) as src:
        # defining new projection.
        dst_crs = 'EPSG:3031'
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        
        # copying the original metadata for the new conversion.
        out_meta = src.meta.copy()
        out_meta.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height,
        })

        # opening a new file for the reprojected data. 
        with rasterio.open(dest_path, 'w', **out_meta) as dest:
            # reprojecting each band.
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dest, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)
                
        # using matplotlib to show the reprojected plot. 
        with rasterio.open(dest_path) as reprojected_src:
            data = reprojected_src.read(1, masked=True)
            fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': ccrs.SouthPolarStereo()})
            ax.coastlines()
            ax.stock_img()
            # defining the extent. 
            extent = (transform[2], transform[2] + transform[0] * width,
                      transform[5] + transform[4] * height, transform[5])
            ax.set_extent(extent, crs=ccrs.SouthPolarStereo())
            im = ax.imshow(data, cmap='viridis', extent=extent,
                           transform=ccrs.SouthPolarStereo(), interpolation='nearest')
            plt.title(f'Reprojected Tiffs (EPSG: 3031) for {year}')
            plt.colorbar(im, ax=ax, orientation='vertical', label='Elevation')
            # saving the plot as png file. 
            plt.savefig(f'Reprojected_Tiffs_{year}.png')
            plt.close()

def fill_gaps(in_path, out_path, year):
    ''' Task 4: filling gaps algorthim, uses rasterio fillnodata to fill the gaps between the flight lines'''
    with rasterio.open(in_path) as src:
        data = src.read(1)
        mask = src.read_masks(1).astype(bool)
        # sets the searching distance for filling gaps. when the smoothing iterations is above 0 the plot turns white so it is set to 0. 
        data_filled = fillnodata(data, mask=mask, max_search_distance=500, smoothing_iterations=0)
        out_meta = src.meta.copy()
        # writing the file to the destination path. 
        with rasterio.open(out_path, "w", **out_meta) as dest:
            # writes to the first band of the destination file.
            dest.write(data_filled, 1)
        # plots using matplotlib using viridis colour scheme. 
        plt.imshow(data_filled, cmap='viridis')
        plt.title(f'Filled Gaps for {year}')
        # creating elevation legend.
        plt.colorbar(label='Elevation')
        # saving figure based on year. 
        plt.savefig(f'Filled_Gaps_{year}.png')
        plt.close()

def parse_arguments():
    ''' This parses arguments instead of hard coding. This is quite complex, It was attempted to combine functions but the code did not run.'''
    parser = argparse.ArgumentParser(description="Process LVIS data.")
    
    subparsers = parser.add_subparsers(dest='command', help='Sub-command help', required=True)

    # Processing files
    process_parser = subparsers.add_parser('process_files', help='Process .h5 files')
    process_parser.add_argument('--input_dir', type=str, required=True, help="Directory with input .h5 files")
    process_parser.add_argument('--scratch_dir', type=str, required=True, help="Directory for temporary outputs")
    process_parser.add_argument('--year', type=str, choices=['2015', '2009'], required=True, help="Specify the year (2015 or 2009)")

    # Merging files
    merge_parser = subparsers.add_parser('merge', help='Merge .tif files')
    merge_parser.add_argument('--scratch_dir', type=str, required=True, help="Directory with temporary .tif files")
    merge_parser.add_argument('--output_tif', type=str, required=True, help="Output merged TIFF file")
    merge_parser.add_argument('--year', type=str, choices=['2015', '2009'], required=True, help="Specify the year (2015 or 2009)")

    # Reprojecting files
    reproject_parser = subparsers.add_parser('reproject', help='Reproject TIFF file')
    reproject_parser.add_argument('--src_path', type=str, required=True, help="Source TIFF file")
    reproject_parser.add_argument('--dest_path', type=str, required=True, help="Destination reprojected TIFF file")
    reproject_parser.add_argument('--year', type=str, choices=['2015', '2009'], required=True, help="Specify the year (2015 or 2009)")
    
    # Filling gaps in TIFFs
    fill_parser = subparsers.add_parser('fill_gaps', help='Fill gaps in TIFF file')
    fill_parser.add_argument('--in_path', type=str, required=True, help="Input TIFF file with gaps")
    fill_parser.add_argument('--out_path', type=str, required=True, help="Output TIFF file with filled gaps")
    fill_parser.add_argument('--year', type=str, choices=['2015', '2009'], required=True, help="Specify the year (2015 or 2009)")
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    
    if args.command == 'process_files':
        process_files(args.input_dir, args.scratch_dir)
    
    # process_files --input_dir /geos/netdata/oosa/assignment/lvis/2015/ --scratch_dir /home/s2478921/scratch/

    elif args.command == 'merge':
        merge_tif(args.scratch_dir, args.output_tif, args.year)

        #python testt3.py merge --scratch_dir /home/s2478921/scratch --output_tif /home/s2478921/2015merged.vrt --year 2015

    elif args.command == 'reproject':
        reproject_to_3031(args.src_path, args.dest_path, args.year)

        #python testt3.py reproject --src_path /home/s2478921/2015merged.vrt --dest_path /home/s2478921/reprojected_output2015.tif --year 2015
    
    elif args.command == 'fill_gaps':
        fill_gaps(args.in_path, args.out_path, args.year)

        #python testt3.py fill_gaps --in_path /home/s2478921/reprojected_output2015.tif --out_path /home/s2478921/filled_output2015.tif --year 2015