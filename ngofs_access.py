# Query NGOF, NEGOF, and NWGOF data over THREDDS

import numpy as np
import pandas as pd
import sys, os
import datetime
import netCDF4
from enum import Enum
from bs4 import BeautifulSoup
import requests
import re



    # Definitions
# NGOFS thredds url string format
ngofs_thredds_directory = "https://opendap.co-ops.nos.noaa.gov/thredds/catalog/NOAA/NGOFS/MODELS/{}{}/catalog.html"
ngofs_nc_prefix = "http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NGOFS/MODELS/"

    # User options
cast = 1 # forecast
targetTime = datetime.datetime(2019, 10, 5, 0, 0, 0)

    # List available NetCDF files for desired day
# Build thredds directory url based on target year, month
url = ngofs_thredds_directory.format(targetTime.year, targetTime.month)
# Request entire directory
page = requests.get(url)
data = page.text
soup = BeautifulSoup(data, features = "lxml")
items = [link.get('href') for link in soup.find_all('a')]
# Filter directory by target day
filterStr = '.{}{}{:02d}.'.format(targetTime.year, targetTime.month, targetTime.day)
items = [i for i in items if filterStr in i]
# Separate into types
itemsNowcast = [i for i in items if "nowcast" in i]
itemsForecast = [i for i in items if "forecast" in i]
itemsFields = [i for i in items if "fields" in i]
# Sort by field (6th) in url that indicates when this file added
itemsForecast.sort(key = lambda x: x.split('.')[6], reverse = True)

    # Open NetCDF connection
# Select most up-to-date forecast for selected day
nc_url = "{}/{}{}/{}".format(ngofs_nc_prefix, str(targetTime.year), str(targetTime.month),
        re.sub(r'.*/', '', itemsForecast[0]))
# Connect to nc data
nc = netCDF4.Dataset(nc_url)

# Print nc overview
print(nc)

# Convert time to nc index
idxTime = netCDF4.date2index(targetTime, nc.variables["time"], select = "nearest")

# Get lat, lon (all)
lat = nc.variables["lat"][:]
lon = nc.variables["lon"][:]
u = nc.variables["u"][idxTime, 0, :]   # Select u at time, surface (0), and all stations (:)
v = nc.variables["v"][idxTime, 0, :]
siglay = nc.variables["siglay"][:]     # Seems to be _where_ sig, idx 0 = -0.01, idx max = -0.98
siglev = nc.variables["siglev"][:]     # Seems to be _what_ sig, where surface is index 0

print(lon.shape)
print(v.shape)

points = pd.DataFrame(list(zip(lat, lon, u, v)), columns = ["lat", "lon", "u", "v"])

print(points)

print(siglay)
print("")
print(siglev)

