
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
import datetime
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

#######################################

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

def convert_longitude(lon):
    return lon + 360 if lon < 0 else lon

def process_files(input_dir, output_dir):
    fileList = glob(os.path.join(input_dir, '*.h5'))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    margin = 2
    mid_x = convert_longitude(-100.000)
    mid_y = -75.167
    shift_right = 2
    shift_down = 0.2
    minX = mid_x - margin + shift_right
    maxX = mid_x + margin + shift_right
    minY = mid_y - margin - shift_down
    maxY = mid_y + margin - shift_down
    tileSize = 2
    xTiles = int((maxX - minX) / tileSize)
    yTiles = int((maxY - minY) / tileSize)
    for index, filePath in enumerate(fileList):
        if index % 3 != 0:
            continue
        b = lvisGround(filePath, onlyBounds=True)
        if b.bounds[0] > maxX or b.bounds[2] < minX or b.bounds[1] > maxY or b.bounds[3] < minY:
            continue
        for i in range(xTiles):
            for j in range(yTiles):
                tileMinX = minX + i * tileSize
                tileMaxX = tileMinX + tileSize
                tileMinY = minY + j * tileSize
                tileMaxY = tileMinY + tileSize
                if b.bounds[0] > tileMaxX or b.bounds[2] < tileMinX or b.bounds[1] > tileMaxY or b.bounds[3] < tileMinY:
                    continue
                lvis = lvisGround(filePath, minX=tileMinX, minY=tileMinY, maxX=tileMaxX, maxY=tileMaxY)
                if lvis.nWaves == 0:
                    continue
                lvis.setElevations()
                lvis.estimateGround()
                base_name = os.path.splitext(os.path.basename(filePath))[0]
                output_filename = os.path.join(output_dir, f'{base_name}_tile_{i}_{j}.tif')
                lvis.writeTiff(data=lvis.zG, x=lvis.lon, y=lvis.lat, res=0.01, filename=output_filename)

def merge_tif(scratch_dir, out_path, year):
    dem_fps = glob(os.path.join(scratch_dir, '*.tif'))
    src_files_to_mosaic = [rasterio.open(fp) for fp in dem_fps]
    mosaic, out_trans = merge(src_files_to_mosaic)
    for src in src_files_to_mosaic:
        src.close()
    plt.imshow(mosaic[0], cmap='terrain')
    plt.savefig(f'Merged_Tiffs_{year}.png')
    plt.close()
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "count": mosaic.shape[0],
        "dtype": mosaic.dtype,
    })
    with rasterio.open(out_path, "w", **out_meta) as dest:
        dest.write(mosaic)

def reproject_to_3031(src_path, dest_path, year):
    with rasterio.open(src_path) as src:
        dst_crs = 'EPSG:3031'
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        out_meta = src.meta.copy()
        out_meta.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height,
        })
        with rasterio.open(dest_path, 'w', **out_meta) as dest:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dest, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)
        with rasterio.open(dest_path) as reprojected_src:
            data = reprojected_src.read(1, masked=True)
            fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': ccrs.SouthPolarStereo()})
            ax.coastlines()
            ax.stock_img()
            extent = (transform[2], transform[2] + transform[0] * width,
                      transform[5] + transform[4] * height, transform[5])
            ax.set_extent(extent, crs=ccrs.SouthPolarStereo())
            im = ax.imshow(data, cmap='terrain', extent=extent,
                           transform=ccrs.SouthPolarStereo(), interpolation='nearest')
            plt.title(f'Reprojected Tiffs (EPSG: 3031) for {year}')
            plt.colorbar(im, ax=ax, orientation='vertical', label='Elevation')
            plt.savefig(f'Reprojected_Tiffs_{year}.png')
            plt.close()

def fill_gaps(in_path, out_path, year):
    with rasterio.open(in_path) as src:
        data = src.read(1)
        mask = src.read_masks(1).astype(bool)
        data_filled = fillnodata(data, mask=mask, max_search_distance=500, smoothing_iterations=0)
        out_meta = src.meta.copy()
        with rasterio.open(out_path, "w", **out_meta) as dest:
            dest.write(data_filled, 1)

        plt.imshow(data_filled, cmap='terrain')
        plt.title(f'Filled Gaps for {year}')
        plt.colorbar(label='Elevation')
        plt.savefig(f'Filled_Gaps_{year}.png')
        plt.close()

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process LVIS data.")
    
    subparsers = parser.add_subparsers(dest='command', help='Sub-command help', required=True)

    # Process files
    process_parser = subparsers.add_parser('process_files', help='Process .h5 files')
    process_parser.add_argument('--input_dir', type=str, required=True, help="Directory with input .h5 files")
    process_parser.add_argument('--scratch_dir', type=str, required=True, help="Directory for temporary outputs")
    process_parser.add_argument('--year', type=str, choices=['2015', '2009'], required=True, help="Specify the year (2015 or 2009)")

    # Merge files
    merge_parser = subparsers.add_parser('merge', help='Merge .tif files')
    merge_parser.add_argument('--scratch_dir', type=str, required=True, help="Directory with temporary .tif files")
    merge_parser.add_argument('--output_tif', type=str, required=True, help="Output merged TIFF file")
    merge_parser.add_argument('--year', type=str, choices=['2015', '2009'], required=True, help="Specify the year (2015 or 2009)")

    # Reproject files
    reproject_parser = subparsers.add_parser('reproject', help='Reproject TIFF file')
    reproject_parser.add_argument('--src_path', type=str, required=True, help="Source TIFF file")
    reproject_parser.add_argument('--dest_path', type=str, required=True, help="Destination reprojected TIFF file")
    reproject_parser.add_argument('--year', type=str, choices=['2015', '2009'], required=True, help="Specify the year (2015 or 2009)")
    
    # Fill gaps in TIFF
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