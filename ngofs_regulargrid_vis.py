from osgeo import gdal, ogr, osr, gdalconst
from gdalconst import GA_ReadOnly
import numpy as np
import cmocean
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def initMap(llx, lly, urx, ury):
    plt.figure(figsize=(12,6))
    m = Basemap(llcrnrlon = llx, llcrnrlat = lly, urcrnrlon = urx, urcrnrlat = ury,
        resolution = "h", epsg = "4269")
    #m.drawmapboundary(fill_color = 'aqua')
    m.fillcontinents(color = '0.8', alpha = 0.5)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,120.,.5), labels = [1,0,0,0], color = 'gray')
    m.drawmeridians(np.arange(-180.,180,1), labels = [0,0,0,1], color = 'gray')

    return m

def dRegion(llx, lly, urx, ury):
    return   { "llx" : llx,
               "lly" : lly,
               "urx" : urx,
               "ury" : ury,
             }


## Options
# Input file directory
indir = "test/data/"
# Input file prefix
prefix = "test"
# Number of bands in rasters (can be less, but not more)
numBands = 8
# Model used to produce forecasts
model = "NWGOFS"
# Select regions of interest. None is full raster region
regions = [None,  dRegion(-95.1, 29.1, -94.4, 29.9) ]

## Map display controls
# Temp constants
titleFmt = "{} forecast ({})\nmodel: {}, time: {}, band: {}\ngenerated: {}"
minWaterSpeed = 0.3
maxWaterSpeed = 2.0
minTemp = 10.0
maxTemp = 20.0
minSalt = 0.0
maxSalt = 36.0

# Setup paths to all rasters
fileFormat = indir + "/" + prefix + "_{}.tiff"
outFormat = indir + "/" + prefix + "_{}_r{}_b{}.png"
file_v = fileFormat.format("v")
file_u = fileFormat.format("u")
file_temp = fileFormat.format("temp")
file_salt = fileFormat.format("salinity")

# Open all rasters
ds_v = gdal.Open(file_v, GA_ReadOnly)
ds_u = gdal.Open(file_u, GA_ReadOnly)
ds_temp = gdal.Open(file_temp, GA_ReadOnly)
ds_salt = gdal.Open(file_salt, GA_ReadOnly)

# Determine lat, lon coordinates of region extent
ulx, xres, xskew, uly, yskew, yres  = ds_v.GetGeoTransform()
lrx = ulx + (ds_v.RasterXSize * xres)
lry = uly + (ds_v.RasterYSize * yres)
llx = lrx - (ds_v.RasterXSize * xres)
lly = lry
urx = ulx + (ds_v.RasterXSize * xres)
ury = uly

r = 0
for region in regions:
    r = r + 1
    if region is None:
        rllx = llx
        rlly = lly
        rurx = urx
        rury = ury
    else:
        rllx = region["llx"]
        rlly = region["lly"]
        rurx = region["urx"]
        rury = region["ury"]

    rulx = rllx
    ruly = rury
    rlrx = rurx
    rlry = rlly

    row1 = int((ruly - uly)/yres)
    col1 = int((rulx - ulx)/xres)
    row2 = int((rlry - uly)/yres)
    col2 = int((rlrx - ulx)/xres)

    for i in range(numBands):
        b = i + 1

        ## Water currents
        plt.close('all')
        # Base map
        m = initMap(rllx, rlly, rurx, rury)
        plt.title(titleFmt.format("Currents", "Knots", model, "???", str(b), "???"))
        # Get data
        data_u = ds_u.GetRasterBand(b).ReadAsArray(col1, row1, col2 - col1 , row2 - row1 )
        data_v = ds_v.GetRasterBand(b).ReadAsArray(col1, row1, col2 - col1 , row2 - row1 )
        data_u[data_u < -99998] = np.nan
        data_v[data_v < -99998] = np.nan
        speed = np.sqrt(data_u * data_u + data_v * data_v)

        x = np.linspace(m.llcrnrx, m.urcrnrx, data_u.shape[1])
        y = np.linspace(m.urcrnry, m.llcrnry, data_u.shape[0])
        xgrid = np.arange(0, data_u.shape[1], 10)
        ygrid = np.arange(0, data_u.shape[0], 10)
        points = np.meshgrid(ygrid, xgrid)
        x = np.linspace(m.llcrnrx, m.urcrnrx, data_u.shape[1])
        y = np.linspace(m.urcrnry, m.llcrnry, data_u.shape[0])
        xx, yy = np.meshgrid(x, y)

        x, y = m(xx, yy)
        contour = m.contourf(x, y, speed, 50, latlon = True, cmap = cmocean.cm.speed,
            vmin = minWaterSpeed, vmax = maxWaterSpeed)
        m.colorbar(contour)
        m.quiver(xx[points], yy[points], data_u[points], data_v[points], alpha = 0.6,
            latlon = True)
        plt.savefig(outFormat.format("currents", str(r), str(b)))

        ## Salinity
        plt.close('all')
        # Base map
        m = initMap(rllx, rlly, rurx, rury)
        plt.title(titleFmt.format("Salinity", "PSU", model, "???", str(b), "???"))
        # Plot
        data = ds_salt.GetRasterBand(b).ReadAsArray(col1, row1, col2 - col1 , row2 - row1 )
        data[data < -99998] = np.nan
        x = np.linspace(m.llcrnrx, m.urcrnrx, data.shape[1])
        y = np.linspace(m.urcrnry, m.llcrnry, data.shape[0])
        xx, yy = np.meshgrid(x, y)
        x, y = m(xx, yy)
        contour = m.contourf(x, y, data, 50, cmap = cmocean.cm.haline,
            vmin = minSalt, vmax = maxSalt, latlon = True)
        m.colorbar(contour)
        plt.savefig(outFormat.format("salinity", str(r), str(b)))

        ## Water temperature
        plt.close('all')
        # Base map
        m = initMap(rllx, rlly, rurx, rury)
        plt.title(titleFmt.format("Water Temperature", "Celsius", model, "???", str(b), "???"))
        # Plot
        data = ds_temp.GetRasterBand(b).ReadAsArray(col1, row1, col2 - col1 , row2 - row1 )
        data[data < -99998] = np.nan
        x = np.linspace(m.llcrnrx, m.urcrnrx, data.shape[1])
        y = np.linspace(m.urcrnry, m.llcrnry, data.shape[0])
        xx, yy = np.meshgrid(x, y)
        x, y = m(xx, yy)
        contour = m.contourf(x, y, data, 50, cmap = cmocean.cm.thermal, vmin = minTemp, vmax = maxTemp, latlon = True)
        m.colorbar(contour)
        plt.savefig(outFormat.format("temp", str(r), str(b)))

