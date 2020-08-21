from pylab import *
import matplotlib.tri as Tri
import netCDF4
import datetime as dt
from optparse import OptionParser
from scipy.interpolate import griddata
from osgeo import gdal, ogr, osr, gdalconst
import rasterio

# Converts a NECOFS NetCDF data product to rasters
# Where each raster holds data from one variable (i.e. water velocity u component)
# and each raster's band corresponds to a discrete time index
# Only stores surface-level data. So, surface currents are captured but not subsurface

def main():

    parser = OptionParser()
    parser.add_option("-n", "--nc",
                      help = "Path to NECOFS NetCDF.",
                      default = "http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc")
    parser.add_option("-s", "--start",
                      help = "Desired time of first forecast. If not set, use current time. Format: Y:M:D:H:M:S",
                      default = None)
    parser.add_option("-v", "--verbose", action = "store_true",
                      help = "Print verbose NetCDF metadata.",
                      default = False)
    parser.add_option("-i", "--info", action = "store_true",
                      help = "Print NetCDF information and exit.",
                      default = False)
    parser.add_option("-b", "--bounds",
                      help = "Comma-delimited spatial bounds for region raster i.e. 'lon_min,lon_max,lat_min,lat_max'",
                      default = "-70.97,-70.82,42.25,42.35")
    parser.add_option("-r", "--rows",
                      help = "Number of rows of target raster",
                      default = 1000)
    parser.add_option("-c", "--cols",
                      help = "Number of columns of target raster",
                      default = 1500)
    parser.add_option("-t", "--times", type = "int",
                      help = "Number of sequential forecast hours where each is a raster band.",
                      default = 3)
    parser.add_option("-g", "--geotiff_prefix",
                      help = "Prefix of file paths to store output rasters",
                      default = None)
    parser.add_option("-p", "--plot",
                      help = "File path to save plot.",
                      default = None)
    (options, args) = parser.parse_args()

    ncUrl = options.nc
    verbose = options.verbose
    infoOnly = options.info
    bounds = [float(b) for b in options.bounds.split(",")]
    outDims = [int(options.rows), int(options.cols)]
    hours = options.times
    startTimeIn = options.start
    outFilePlot = options.plot
    outFileRasterPre = options.geotiff_prefix

    print(ncUrl)

    nc = netCDF4.Dataset(ncUrl)
    ncVars = nc.variables

    print("NetCDF: {t}".format(t = nc.title))
    print("------------------")
    print("URL: {url}".format(url = ncUrl))
    ncStartTime = netCDF4.num2date(ncVars["time"][0], ncVars["time"].units)
    ncEndTime = netCDF4.num2date(ncVars["time"][-1], ncVars["time"].units)
    print("Duration: {s} ---- {e}".format(s = ncStartTime, e = ncEndTime))

    if verbose == False:
        print("Keys:")
        print(ncVars.keys())
    else:
        print("Keys:")
        print("----------")
        for v in ncVars.keys():
            print("Key: {v}".format(v = v))
            print(ncVars[v])
            print("")
        print(nc)

    if infoOnly:
        exit(0)

    # Output raster info
    print("\nOutput rasters:")
    print("Rows: {r}, columns: {c}, bands: {b}".format(r = outDims[0], c = outDims[1], b = hours))

    # Desired access time

    if startTimeIn is not None:
        startTime = dt.datetime.strptime(startTimeIn, '%Y:%m:%d:%H:%M:%S')
    else:
        startTime = dt.datetime.utcnow() + dt.timedelta(hours=10)
    itime = netCDF4.date2index(startTime, ncVars["time"], select = "nearest")
    dtime = netCDF4.num2date(ncVars["time"][itime], ncVars["time"].units)
    daystr = dtime.strftime('%Y-%b-%d %H:%M')
    print("Desired start time: {s})".format(s = startTime))
    print("Desired number of hourly bands: {b}".format(b= hours))
    print("Forecast band times:")
    idxs = list(range(itime, itime + hours))
    for i in range(len(idxs)):
        dtime = netCDF4.num2date(ncVars["time"][idxs[i]], ncVars["time"].units)
        daystr = dtime.strftime('%Y-%b-%d %H:%M')
        print("    [{i}] {s}".format(i = i, s = daystr))

    # Get lon,lat coordinates for nodes (depth)
    lat = nc['lat'][:]
    lon = nc['lon'][:]

    # Get lan, lon coords for cell centers
    latc = nc['latc'][:]
    lonc = nc['lonc'][:]

    # Only select points within bounds
    ax = bounds
    ind = argwhere((lon >= ax[0]) & (lon <= ax[1]) & (lat >= ax[2]) & (lat <= ax[3]))
    indc = argwhere((lonc >= ax[0]) & (lonc <= ax[1]) & (latc >= ax[2]) & (latc <= ax[3]))
    lat = lat[ind][:, 0]
    lon = lon[ind][:, 0]
    latc = latc[indc][:, 0]
    lonc = lonc[indc][:, 0]

    # Init meshgrid for interpolation
    xx, yy = np.mgrid[bounds[0]:bounds[1]:complex(0, outDims[1]), bounds[3]:bounds[2]:complex(0, outDims[0])]

    # Init plot
    fig, axs = plt.subplots(len(idxs), 2, figsize=(20, 6 * len(idxs)))
    axs[0][1].set_title("Salinity (ppt)", y = 1.08, fontsize = 10)
    axs[0][0].set_title("Depth (meters)", y = 1.08, fontsize = 10)

    # Init output arrays
    hout = np.zeros((outDims[0], outDims[1]))
    uout = np.zeros((outDims[0], outDims[1], len(idxs)))
    vout = np.zeros((outDims[0], outDims[1], len(idxs)))
    sout = np.zeros((outDims[0], outDims[1], len(idxs)))

    # Get data
    h = ncVars["h"][:]
    h = h[ind][:, 0]
    hgrid = griddata((lon, lat), h, (xx, yy), method = "cubic").T

    hout = hgrid.copy()

    i = 0
    for itime in idxs:
        dtime = netCDF4.num2date(ncVars["time"][itime], ncVars["time"].units)
        daystr = dtime.strftime('%Y-%b-%d %H:%M')

        u = ncVars["u"][itime, 0, :]
        u = u[indc][:, 0]
        ugrid = griddata((lonc, latc), u, (xx, yy), method = "cubic").T
        uout[:, :, i] = ugrid.copy()

        v = ncVars["v"][itime, 0, :]
        v = v[indc][:, 0]
        vgrid = griddata((lonc, latc), v, (xx, yy), method = "cubic").T
        vout[:, :, i] = vgrid.copy()

        s = ncVars["salinity"][itime, 0, :]
        s = s[ind][:, 0]
        sgrid = griddata((lon, lat), s, (xx, yy), method = "cubic").T
        sout[:, :, i] = sgrid.copy()

        im = axs[i][0].imshow(hgrid, cmap = "terrain_r")
        plt.colorbar(im, ax = axs[i][0])
        skip = int(pow(outDims[0] * outDims[0] + outDims[1] * outDims[1], 0.5 ) * 0.015)
        Q = axs[i][0].quiver(range(0, outDims[1], skip), range(0, outDims[0], skip), ugrid[::skip, ::skip], vgrid[::skip, ::skip], scale = 20)
        keyVelocity = 0.5
        keyStr = 'Water: %3.1f m/s' % keyVelocity
        qk = axs[i][0].quiverkey(Q, 1.05, 1.06, keyVelocity, keyStr, labelpos = 'W')

        im = axs[i][1].imshow(sgrid, cmap = "bone")
        plt.colorbar(im, ax = axs[i][1])
        fig.suptitle("Model: {}, Start: {}".format("NECOFS", daystr), fontsize = 12)

        xtickPos = list(range(0, outDims[1], int(outDims[1] / 4)))
        xtickLabel = ["{:.2f}°W".format(s) for s in linspace(bounds[0], bounds[1], 4)]
        ytickPos = list(range(0, outDims[0], int(outDims[0] / 4)))
        ytickLabel = ["{:.2f}°N".format(s) for s in linspace(bounds[3], bounds[2], 4)]
        axs[i][0].set_yticks(ytickPos)
        axs[i][0].set_yticklabels(ytickLabel, fontsize = 8)
        axs[i][1].set_yticks(ytickPos)
        axs[i][1].set_yticklabels(ytickLabel, fontsize = 8)
        axs[i][0].set_xticks(xtickPos)
        axs[i][0].set_xticklabels(xtickLabel, fontsize = 8)
        axs[i][1].set_xticks(xtickPos)
        axs[i][1].set_xticklabels(xtickLabel, fontsize = 8)

        i += 1

    # Plot data
    if outFilePlot is not None:
        plt.savefig(outFilePlot)
    else:
        plt.show()

    if outFileRasterPre is not None:
        # Store raster
        hOutFileRaster = outFileRasterPre + "_height.tiff"
        uOutFileRaster = outFileRasterPre + "_uwater.tiff"
        vOutFileRaster = outFileRasterPre + "_vwater.tiff"
        sOutFileRaster = outFileRasterPre + "_salinity.tiff"

        # CRS source: http://easterndivision.s3.amazonaws.com/Marine/MooreGrant/5_SMAST_SASI_Oceanograhic_variability.pdf
        crs = "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
        transform = rasterio.transform.from_bounds(bounds[0], bounds[2], bounds[1], bounds[3], outDims[1], outDims[0])

        h_ds = rasterio.open(hOutFileRaster, 'w', driver = "GTiff", height = outDims[0], width = outDims[1],
                             count = 1, dtype = str(hout.dtype), crs = crs, transform = transform)
        h_ds.write(hout, 1)
        h_ds.close()

        u_ds = rasterio.open(uOutFileRaster, 'w', driver = "GTiff", height = outDims[0], width = outDims[1],
                             count = len(idxs), dtype = str(uout.dtype), crs = crs, transform = transform)
        u_ds.write(np.rollaxis(uout, axis = 2))
        u_ds.close()

        v_ds = rasterio.open(vOutFileRaster, 'w', driver = "GTiff", height = outDims[0], width = outDims[1],
                             count = len(idxs), dtype = str(vout.dtype), crs = crs, transform = transform)
        v_ds.write(np.rollaxis(vout, axis = 2))
        v_ds.close()

        s_ds = rasterio.open(sOutFileRaster, 'w', driver = "GTiff", height = outDims[0], width = outDims[1],
                             count = len(idxs), dtype = str(sout.dtype), crs = crs, transform = transform)
        s_ds.write(np.rollaxis(sout, axis = 2))
        s_ds.close()

if __name__ == '__main__':
    main()

