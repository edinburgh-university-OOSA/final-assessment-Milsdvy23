import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os

# Define directories
directory_2009 = '/geos/netdata/oosa/assignment/lvis/2009/'
directory_2015 = '/geos/netdata/oosa/assignment/lvis/2015/'
output_directory = '/home/s2478921/oose/final-assessment-Milsdvy23/hdf'

# Files to exclude
exclude_files = set([
    'ILVIS1B_AQ2009_1020_R1408_049700.h5',
    'ILVIS1B_AQ2009_1020_R1408_052195.h5',
    'ILVIS1B_AQ2009_1020_R1408_053614.h5',
    'ILVIS1B_AQ2009_1020_R1408_068453.h5',
    'ILVIS1B_AQ2009_1020_R1408_070488.h5',
    'ILVIS1B_AQ2009_1020_R1408_071909.h5',
    'ILVIS1B_AQ2015_1012_R1605_048197.h5',
    'ILVIS1B_AQ2015_1012_R1605_055228.h5',
    'ILVIS1B_AQ2015_1012_R1605_070498.h5',
    'ILVIS1B_AQ2015_1017_R1605_043439.h5',
    'ILVIS1B_AQ2015_1017_R1605_055462.h5',
    'ILVIS1B_AQ2015_1017_R1605_067952.h5',
    'ILVIS1B_AQ2015_1017_R1605_069264.h5',
    'ILVIS1B_AQ2015_1017_R1605_071670.h5',
    'ILVIS1B_AQ2015_1017_R1605_082081.h5',
    'ILVIS1B_AQ2015_1027_R1605_048105.h5',
    'ILVIS1B_AQ2015_1027_R1605_051092.h5',
    'ILVIS1B_AQ2015_1027_R1605_058362.h5',
    'ILVIS1B_AQ2015_1027_R1605_059491.h5',
    'ILVIS1B_AQ2015_1027_R1605_060191.h5',
    'ILVIS1B_AQ2015_1027_R1605_072553.h5',
    'ILVIS1B_AQ2015_1027_R1605_073584.h5',
    'ILVIS1B_AQ2015_1027_R1605_074433.h5',
    'ILVIS1B_AQ2015_1027_R1605_075604.h5',
    'ILVIS1B_AQ2015_1027_R1605_079057.h5'
])

# Filter files
files_2009 = [
    f for f in os.listdir(directory_2009) if f.endswith('.h5') and f not in exclude_files
]

files_2015 = [
    f for f in os.listdir(directory_2015) if f.endswith('.h5') and f not in exclude_files
]

def extract_data(directory, files, lat_key, lon_key):
    latitudes, longitudes = [], []
    for filename in files:
        filepath = os.path.join(directory, filename)
        with h5py.File(filepath, 'r') as f:
            lat_start = f['LAT0'][:]
            lat_end = f[lat_key][:]
            lon_start = f['LON0'][:]
            lon_end = f[lon_key][:]
            
            midpoint_lat = (lat_start + lat_end) / 2
            midpoint_lon = (lon_start + lon_end) / 2

            # Sample every other point
            latitudes.extend(midpoint_lat[::2])
            longitudes.extend(midpoint_lon[::2])
    
    return latitudes, longitudes

# Combine data
latitudes_2009, longitudes_2009 = extract_data(directory_2009, files_2009, 'LAT1023', 'LON1023')
latitudes_2015, longitudes_2015 = extract_data(directory_2015, files_2015, 'LAT527', 'LON527')

# Plot combined map
plt.figure(figsize=(12, 8))
m = Basemap(projection='cyl',
            llcrnrlat=min(latitudes_2009 + latitudes_2015) - 5,
            urcrnrlat=max(latitudes_2009 + latitudes_2015) + 5,
            llcrnrlon=min(longitudes_2009 + longitudes_2015) - 5,
            urcrnrlon=max(longitudes_2009 + longitudes_2015) + 5,
            resolution='l')
m.drawcoastlines()
m.drawcountries()
m.drawparallels(range(-90, 91, 10), labels=[1, 0, 0, 0])
m.drawmeridians(range(-180, 181, 20), labels=[0, 0, 0, 1])

x_2009, y_2009 = m(longitudes_2009, latitudes_2009)
m.scatter(x_2009, y_2009, s=1, c='red', label='2009 Data')

x_2015, y_2015 = m(longitudes_2015, latitudes_2015)
m.scatter(x_2015, y_2015, s=1, c='blue', label='2015 Data')

plt.title('Combined Map of Remaining 2009 and 2015 HDF5 Files')
plt.legend()
plt.savefig(os.path.join(output_directory, "combined_map.png"))
plt.show()
