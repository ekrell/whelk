# Query NGOF, NEGOF, and NWGOF data over THREDDS)

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

def hour2interval(hour, intervals):
    # Returns the update time (in hours) to select closest previous product
    # If None is returned, then the previous day should be used
    interval = None
    if hour < intervals[0]:
        return interval
    for i in intervals:
        if i <= hour:
            interval = i
        else:
            break
    return interval

    # Definitions
# Models (NGOFS and its nested models)
models = [ { "name"      : "ngofs",
             "parent"    : None,
             "thredds"   : "https://opendap.co-ops.nos.noaa.gov/thredds/catalog/NOAA/NGOFS/MODELS/",
             "ncPrefix" : "http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NGOFS/MODELS/",
           },
           { "name"      : "nwgofs",
             "parent"    : "ngofs",
             "thredds"   : "https://opendap.co-ops.nos.noaa.gov/thredds/catalog/NOAA/NWGOFS/MODELS/",
             "ncPrefix" : "http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NWGOFS/MODELS/",

           },
           { "name"      : "negofs",
             "parent"    : "ngofs",
             "thredds"   : "https://opendap.co-ops.nos.noaa.gov/thredds/catalog/NOAA/NEGOFS/MODELS/",
             "ncPrefix" : "http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NEGOFS/MODELS/",
           },
]
# Time intervals between NGOFS forecasts
intervals = [3, 9, 15, 21]
# Data product types
products = ["nowcast", "forecast", "regulargrid", "field"]
# NC file format
ncFileFmt = "nos.{}.stations.{}.{:04}{:02}{:02d}.t{:02d}z.nc"

    # User options
targetModelIdx = 1;
targetProduct = products[1]
# Time of data product
targetStartTime = datetime.datetime(2019, 10, 2, 20, 0, 0)
targetStopTime =  datetime.datetime(2019, 10, 3, 12, 0, 0)
deltaTime_s = 60 * 10
writeDir = "."

print("--NetCDF Forecast Access--")

# Ensure valid times
if targetStopTime < targetStartTime:
    print("Start time is later than end time!")
    exit(1)
# Build target times
targetTimes = [targetStartTime]
while targetTimes[-1] < targetStopTime:
    targetTimes.append(targetTimes[-1] + datetime.timedelta(seconds = deltaTime_s))
targetTimes.append(targetStopTime)

print("Current time:")
print("    System time: {}".format(datetime.datetime.now()))
print("Requested time duration:")
print("    Start: {}".format(targetStartTime))
print("    Stop: {}".format(targetStopTime))
print("    Delta t: {}".format(datetime.timedelta(seconds = deltaTime_s)))
print("    Duration: {}".format(targetStopTime - targetStartTime))
print("    Num discrete intervals: {}".format(len(targetTimes)))

    # Build data product request
# Select data product available for requested time
# (Closest preceeding, since we can simulate that the future not yet available)
interval = hour2interval(targetStartTime.hour, intervals)
# Handle case where the target hour is before any product have been made for that day
productStartTime = targetStartTime
if interval is None:
    productStartTime = productStartTime - datetime.timedelta(1)
    interval = max(intervals)
# Build url to NetCDF file
ncFile = ncFileFmt.format(models[targetModelIdx]["name"], targetProduct, productStartTime.year, productStartTime.month, productStartTime.day, interval)
ncUrl = "{}/{:02d}{:02d}/{}".format(models[targetModelIdx]["ncPrefix"], productStartTime.year, productStartTime.month, ncFile)

# Output table file
outFileTable = "{}__{}_{}.csv".format(ncFile.replace(".nc", ""), targetStartTime.strftime("%Y%m%d-%H%M%S"), targetStopTime.strftime("%Y%m%d-%H%M%S"))

print("NetCDF product request:")
print("    Start: {}".format(productStartTime))
print("    Forecast interval: {}".format(interval))
print("    Forecast url: {}".format(ncUrl))
print("    Write table as: {}/{}".format(writeDir, outFileTable))

    # Open NetCDF connection
# Connect to nc data
nc = netCDF4.Dataset(ncUrl)
# Check dataset time
ncStartTime = netCDF4.num2date(nc.variables["time"][0], nc.variables["time"].units)
ncEndTime = netCDF4.num2date(nc.variables["time"][-1], nc.variables["time"].units)

print("NetCDF product summary:")
print(nc)
print("------------")
print("NetCDF product time")
print("    Time units: {}".format(nc.variables["time"].units))
print("    Start: {}".format(ncStartTime))
print("    Stop: {}".format(ncEndTime))
print("    Duration: {}".format(ncEndTime - ncStartTime))
print("    Num discrete intervals: {}".format(len(nc.variables["time"][:])))

print("Time sanity check:")
print("    Requested start >= NetCDF start: {}".format(targetStartTime - ncStartTime >= datetime.timedelta(seconds = 0)))
print("    Requested stop  <= NetCDF stop: {}".format(ncEndTime - targetStopTime >= datetime.timedelta(seconds = 0)))

# Convert requested time intervals to nc indices
idxTimes = netCDF4.date2index(targetTimes, nc.variables["time"], select = "nearest")

# Spatial (non-temporal)
name_station = nc.variables["name_station"][:]
names = ["".join([c.decode('utf-8') for c in name])[0:-1] for name in name_station]
lat = nc.variables["lat"][:]
lon = nc.variables["lon"][:]
x = nc.variables["x"][:]
y = nc.variables["y"][:]
h = nc.variables["h"][:]

    # Initialize table (pandas)
# Detetmine number of rows
numStations = len(name_station)
numTimeIntervals = len(idxTimes)
numRows = len(name_station) * len(idxTimes)
df = pd.DataFrame(index = range(numRows), columns = ["station", "lat", "lon", "y", "x", "h", "u", "v", "uwater", "vwater", "temp", "salinity", "zeta"])

points = []
for idxTime in idxTimes:
    timeSlice = netCDF4.num2date(nc.variables["time"][idxTime], nc.variables["time"].units)

    # Temporal, depth
    u = nc.variables["u"][idxTime, 0, :]   # Select u at time, surface (0), and all stations (:)
    v = nc.variables["v"][idxTime, 0, :]
    temp = nc.variables["temp"][idxTime, 0, :]
    salinity = nc.variables["salinity"][idxTime, 0, :]
    # Temporal, no depth
    uwind_speed = nc.variables["uwind_speed"][idxTime, :]
    vwind_speed = nc.variables["vwind_speed"][idxTime, :]
    zeta = nc.variables["zeta"][idxTime, :]

    points.append(pd.DataFrame(list(zip(names, [timeSlice for i in range(numStations)], lat, lon, y, x, h, u, v, uwind_speed, vwind_speed, temp, salinity, zeta)),
        columns = ["station", "time", "lat", "lon", "y", "x", "bathy", "u_water", "v_water", "u_wind", "v_wind", "temp", "salinity", "surface_elev"]))

df = pd.concat(points)

# Write table
df.to_csv("{}/{}".format(writeDir, outFileTable))

