#!/usr/bin/python3

# Query NEGOF, NWGOF regulargrid data over THREDDS

from math import modf
import numpy as np
import pandas as pd
import sys, os
import datetime
import netCDF4
from enum import Enum
from bs4 import BeautifulSoup
import requests
import re
from osgeo import gdal
import osr

models = [ {
        "name"      : "nwgofs",
        "thredds"   : "https://opendap.co-ops.nos.noaa.gov/thredds/catalog/NOAA/NWGOFS/MODELS/",
        "ncPrefix" : "http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NWGOFS/MODELS/",
    },
    {
        "name"      : "negofs",
        "thredds"   : "https://opendap.co-ops.nos.noaa.gov/thredds/catalog/NOAA/NEGOFS/MODELS/",
        "ncPrefix" : "http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NEGOFS/MODELS/",
    },
]

ncFileFmt = "nos.{}.{}.{:04}{:02}{:02d}.t{:02d}z.nc"

rasterFile_u = "test_u.tiff"
rasterFile_v = "test_v.tiff"
rasterFile_temp = "test_temp.tiff"
rasterFile_salt = "test_salt.tiff"

modelIdx = 0
productType = "regulargrid"
productLabel = "f"
year = 2020
month = 1
day = 28
hour = 21    # (3 | 9 | 15 | 21 )
productTime = datetime.datetime(year, month, day, hour, 0, 0)

# Generate list of needed forecast grids, with their netCDF access information
forecastRange = range(0, 8)  # Max: 48

forecasts = [ dict({"n" : None, "ncFile" : None, "ncUrl" : None}) for f in forecastRange ]
for i in range(len(forecastRange)):
    forecasts[i]["n"] = forecastRange[i]
    product = "{}.{}{:03d}".format(productType, productLabel, forecasts[i]["n"])
    forecasts[i]["ncFile"] = ncFileFmt.format(models[modelIdx]["name"], product, productTime.year, productTime.month,
            productTime.day, productTime.hour)
    forecasts[i]["ncUrl"] = "{}/{:02d}{:02d}/{}".format(models[modelIdx]["ncPrefix"], productTime.year,
            productTime.month, forecasts[i]["ncFile"])
    print(forecasts[i])

# Access first forecast
nc = netCDF4.Dataset(forecasts[0]["ncUrl"])

# Get coordinates grids
lat = nc.variables["Latitude"][:]
lon = nc.variables["Longitude"][:]

# Check raster details
(ny, nx) = lat.shape
numBands = len(forecasts)
xmin, ymin, xmax, ymax = [lon.min(), lat.min(), lon.max(), lat.max()]
xres = (xmax - xmin) / float(nx)
yres = (ymax - ymin) / float(ny)
geotransform = (xmin, xres, 0, ymax, 0, -yres)

print("")
print("Extent & resolution:")
print("xmin: {}, xmax: {}, xres: {}, ymin: {}, ymax: {}, yres: {}".format(
    xmin, xmax, xres, ymin, ymax, yres))
print("Coordinate projection: {}".format(nc.CoordinateProjection))
print("")


for name in nc.ncattrs():
    print(name)

# Generate csv header
header = ""
header = header + "# Title: {}\n".format(nc.title)
header = header + "# NC Start: {}\n".format(forecasts[0]["ncFile"])
header = header + "# Institution: {}\n".format(nc.institution)
header = header + "# Source: {}\n".format(nc.source)
header = header + "# History: {}\n".format(nc.history)
header = header + "# References: {}\n".format(nc.references)
header = header + "# CoordinateSystem: {}\n".format(nc.CoordinateSystem)
header = header + "# CoordinateProjection: {}\n".format(nc.CoordinateProjection)
print(header)

# Write csv header

exit()

# Setup NAD83 projection
srs = osr.SpatialReference()
srs.ImportFromEPSG(4269) #NAD83

# Initialize rasters
drv = gdal.GetDriverByName("GTiff")
# Eastward sea water velocity (meters/second)
ds_u = drv.Create(rasterFile_u, nx, ny, numBands, gdal.GDT_Float32)
ds_u.SetGeoTransform(geotransform)
ds_u.SetProjection(srs.ExportToWkt())

# Northward sea water velocity (meters/second)
ds_v = drv.Create(rasterFile_v, nx, ny, numBands, gdal.GDT_Float32)
ds_v.SetGeoTransform(geotransform)
ds_v.SetProjection(srs.ExportToWkt())

# Potential temperature (celsius)
ds_temp = drv.Create(rasterFile_temp, nx, ny, numBands, gdal.GDT_Float32)
ds_temp.SetGeoTransform(geotransform)
ds_temp.SetProjection(srs.ExportToWkt())

# Salinity (PSU)
ds_salt = drv.Create(rasterFile_salt, nx, ny, numBands, gdal.GDT_Float32)
ds_salt.SetGeoTransform(geotransform)
ds_salt.SetProjection(srs.ExportToWkt())


for f in range(len(forecasts)):
    # Access forecast
    nc = netCDF4.Dataset(forecasts[f]["ncUrl"])

    # Band index
    b = f + 1


#
#    # Eastward sea water velocity
#    grid = nc.variables["u_eastward"][:][0][0]
#    band = ds_u.GetRasterBand(b)
#    band.SetNoDataValue(-99999.0)
#    band.WriteArray(grid[::-1,...])
#
#    # Northward sea water velocity
#    grid = nc.variables["v_northward"][:][0][0]
#    band = ds_v.GetRasterBand(b)
#    band.SetNoDataValue(-99999.0)
#    band.WriteArray(grid[::-1,...])
#
#    # Temperature
#    grid = nc.variables["temp"][:][0][0]
#    band = ds_temp.GetRasterBand(b)
#    band.SetNoDataValue(-99999.0)
#    band.WriteArray(grid[::-1,...])
#
#    # Salinity
#    grid = nc.variables["salt"][:][0][0]
#    band = ds_salt.GetRasterBand(b)
#    band.SetNoDataValue(-99999.0)
#    band.WriteArray(grid[::-1,...])

# Write rasters
ds_u.FlushCache()
ds_v.FlushCache()
ds_temp.FlushCache()
ds_salt.FlushCache()

