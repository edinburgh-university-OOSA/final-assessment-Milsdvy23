from osgeo import gdal, gdalconst
import numpy as np
import argparse
import rasterio
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap

def parse_arguments():
    ''' parsing arguments function '''
    parser = argparse.ArgumentParser(description="Calculate volume change between two DEMs and create a GeoTIFF of elevation change.")
    parser.add_argument('input_file1', type=str, help='Path to the first GeoTIFF file')
    parser.add_argument('input_file2', type=str, help='Path to the second GeoTIFF file')
    parser.add_argument('output_file', type=str, help='Output path for the elevation change GeoTIFF')
    return parser.parse_args()

def read_tiff(filename):
    '''reading the tiff to ram depending on filename provided'''
    ds = gdal.Open(filename)
    if not ds:
        raise IOError(f"Unable to open file: {filename}")
    # reading the first band of data into a numpy array.
    data = ds.GetRasterBand(1).ReadAsArray()
    transform = ds.GetGeoTransform()
    # extracting the spatial information. 
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    return data, xOrigin, yOrigin, pixelWidth, pixelHeight, ds.RasterXSize, ds.RasterYSize

def resample_to_match(ds, match_ds):
    ''' Resampling the second geotiff so the resolution and projection matches the first'''
    # creates temporary resampled dataset to match with the new ds dataset,
    resampled = gdal.GetDriverByName('MEM').Create(
        # uses GDT_Float32 to work with elevation data.
        '', match_ds.RasterXSize, match_ds.RasterYSize, 1, gdalconst.GDT_Float32)
    # sets the spatial data to be the same as match_ds
    resampled.SetGeoTransform(match_ds.GetGeoTransform())
    resampled.SetProjection(match_ds.GetProjection())
    # Using gdal.ReprojectImage to reproject and resample ds into the resampled dataset.
    gdal.ReprojectImage(ds, resampled, ds.GetProjection(), match_ds.GetProjection(), gdalconst.GRA_Bilinear)
    # return resampled dataset that will now be viable for calcuating overap. 
    return resampled

def calculate_overlap(data1, data2, metadata1, metadata2):
    ''' Extracts overlapping regions from the two tiffs, using metadata to plot them. '''
    # draws from metadata of both datasets. 
    _, xOrigin1, yOrigin1, pixelWidth1, pixelHeight1, _, _ = metadata1
    _, xOrigin2, yOrigin2, pixelWidth2, pixelHeight2, _, _ = metadata2
    # calcuating overlap bounds,
    x_start = max(xOrigin1, xOrigin2)
    y_start = min(yOrigin1, yOrigin2)
    x_end = min(xOrigin1 + data1.shape[1] * pixelWidth1, xOrigin2 + data2.shape[1] * pixelWidth2)
    y_end = max(yOrigin1 + data1.shape[0] * pixelHeight1, yOrigin2 + data2.shape[0] * pixelHeight2)
    # slicing data arrays by calcuating offset indices. 
    x_offset1 = int((x_start - xOrigin1) / pixelWidth1)
    y_offset1 = int((yOrigin1 - y_start) / abs(pixelHeight1))
    x_offset2 = int((x_start - xOrigin2) / pixelWidth2)
    y_offset2 = int((yOrigin2 - y_start) / abs(pixelHeight2))
    # working out overlap sizes. 
    nX_overlap = int((x_end - x_start) / pixelWidth1)
    nY_overlap = int((y_start - y_end) / abs(pixelHeight1))
    # getting overlap regions. 
    overlap_data1 = data1[y_offset1:y_offset1 + nY_overlap, x_offset1:x_offset1 + nX_overlap]
    overlap_data2 = data2[y_offset2:y_offset2 + nY_overlap, x_offset2:x_offset2 + nX_overlap]
    # returning overlap data and pixel size.
    return overlap_data1, overlap_data2, pixelWidth1, pixelHeight1

def calculate_volume_change(data1, data2, pixelWidth, pixelHeight):
    ''' calcuating volume change by comparing the 2009 and 2015 datasets'''
    # calcuating volume change by minusing the 2015 data from the 2009 data. 
    difference = data2 - data1
    area_per_pixel = abs(pixelWidth * pixelHeight)
    # calcuating the volume change by summing the elevation difference and mutliping by the pixel area. 
    volume_change = np.sum(difference) * area_per_pixel
    return volume_change

def write_tiff(output_filename, data, input_ds):
    ''' creates a geotiff by using the resampled tiff to create a volume change geotiff. '''
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_filename, input_ds.RasterXSize, input_ds.RasterYSize, 1, gdalconst.GDT_Float32)
    out_ds.SetGeoTransform(input_ds.GetGeoTransform())
    out_ds.SetProjection(input_ds.GetProjection())
    out_ds.GetRasterBand(1).WriteArray(data)
    out_ds.FlushCache()

def save_as_png(input_tiff, output_png):
    with rasterio.open(input_tiff) as src:
        data = src.read(1)  # Only reads the fist band. 
        min_val, max_val = data.min(), data.max()

        # Normalising data for color mapping
        norm = Normalize(vmin=min_val, vmax=max_val)
        cmap = get_cmap('viridis')

        # Plotting data using matplotlib. 
        plt.figure(figsize=(10, 8))
        plt.imshow(data, cmap=cmap, norm=norm)
        plt.colorbar(label='Elevation Change')
        plt.title('Elevation Change from 2009 to 2015')
        plt.axis('off')
        plt.savefig(output_png, bbox_inches='tight', pad_inches=0.1)
        plt.close()

def main():
    ''' creating the final output.'''
    args = parse_arguments()

    # Read and processing the first TIFF
    ds1 = gdal.Open(args.input_file1)
    data1, xOrigin1, yOrigin1, pixelWidth1, pixelHeight1, nX1, nY1 = read_tiff(args.input_file1)
    metadata1 = (data1, xOrigin1, yOrigin1, pixelWidth1, pixelHeight1, nX1, nY1)

    # Reading and resampling the second TIFF to match the first.
    ds2 = gdal.Open(args.input_file2)
    resampled_ds2 = resample_to_match(ds2, ds1)
    data2 = resampled_ds2.GetRasterBand(1).ReadAsArray()
    metadata2 = (data2, xOrigin1, yOrigin1, pixelWidth1, pixelHeight1, nX1, nY1)

    # Calculating overlap and volume change
    overlap_data1, overlap_data2, pixelWidth, pixelHeight = calculate_overlap(data1, data2, metadata1, metadata2)
    if overlap_data1.shape == overlap_data2.shape:
        volume_change = calculate_volume_change(overlap_data1, overlap_data2, pixelWidth, pixelHeight)
        print(f"Total Ice Volume Change: {volume_change} cubic meters")
    else:
        print("Error: Unable to match the overlap region shapes.")

    # Creating elevation difference TIFF
    elevation_difference = data2 - data1
    write_tiff(args.output_file, elevation_difference, ds1)
    print(f"Elevation change GeoTIFF saved as: {args.output_file}")

    # Saving as png
    png_output_file = args.output_file.replace('.tif', '.png')
    save_as_png(args.output_file, png_output_file)
    print(f"Elevation change PNG saved as: {png_output_file}")

if __name__ == "__main__":
    main()
