from datetime import datetime, timedelta
import datetime
import time
import sys
import os, shutil
from ftplib import FTP
import logging
import numpy as np
from itertools import groupby
import tempfile, shutil, sys
import calendar
from netCDF4 import Dataset


def downloadFile(grace_data_dir):

    ftp = FTP("podaac.jpl.nasa.gov")
    ftp.login()
    ftp.cwd("/allData/tellus/L3/mascon/RL05/JPL/CRI/netcdf/")

    files = ftp.nlst()
    now_str = str(datetime.datetime.now())
    file_name = now_str+" grace.nc"

    if not os.path.exists(grace_data_dir):
        os.makedirs(grace_data_dir)

    grace_file_path = os.path.join(grace_data_dir, file_name)

    for filename in files:
        if filename.startswith("GRCTellus") and filename.endswith(".nc"):
            local_filename = os.path.join(grace_data_dir, file_name)
            local_file = open(local_filename, 'wb')
            print "Downloading GRACE file..."
            ftp.retrbinary('RETR ' + filename, local_file.write)

            local_file.close()

    return grace_file_path

#grace_data_dir = "/home/tethys/GRACE/global/"
#grace_file_path = downloadFile(grace_data_dir)