import fiona
from netCDF4 import Dataset
import os
from datetime import datetime,timedelta
import calendar
import numpy as np
import tempfile, shutil,sys
import gdal
import ogr
import osr
import requests
import json
import functools
import shapely.geometry #Need this to find the bounds of a given geometry
import shapely.ops
import geojson
import pyproj
from pyproj import Proj, transform
from config import GLOBAL_NETCDF_DIR

#Check if the user is superuser or staff. Only the superuser or staff have the permission to add and manage watersheds.
def user_permission_test(user):
    return user.is_superuser or user.is_staff

def get_global_nc():
    ncfile = None
    for file in os.listdir(GLOBAL_NETCDF_DIR):
        if file.startswith('GRC') and file.endswith('.nc'):
            ncfile = os.path.join(GLOBAL_NETCDF_DIR,file)

    return ncfile

def create_global_tiff(file_name,output_dir,var_name):

    output_dir = os.path.join(output_dir, '')

    xsize, ysize, GeoT, NDV = get_netcdf_info_global(file_name,var_name)


    start_date = '01/01/2002'

    nc_fid = Dataset(file_name, 'r')  # Reading the netcdf file
    nc_var = nc_fid.variables #Get the netCDF variables
    nc_var.keys() #Getting variable keys

    time = nc_var['time'][:]

    date_str = datetime.strptime(start_date, "%m/%d/%Y")  # Start Date string.

    for timestep, v in enumerate(time):

        current_time_step = nc_var[var_name][timestep, :, :]  # Getting the index of the current timestep

        end_date = date_str + timedelta(days=float(v))  # Actual human readable date of the timestep

        ts_file_name = end_date.strftime("%Y_%m_%d")  # Changing the date string format

        data = nc_var[var_name][timestep,:,:]
        data = data[::-1, :]
        driver = gdal.GetDriverByName('GTiff')
        DataSet = driver.Create(output_dir+ts_file_name+'.tif',xsize,ysize,1, gdal.GDT_Float32)
        DataSet.SetGeoTransform([0.0, 0.5, 0.0, 90.0, 0.0, -0.5])
        srs=osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        DataSet.SetProjection(srs.ExportToWkt())

        DataSet.GetRasterBand(1).WriteArray(data)
        DataSet.GetRasterBand(1).SetNoDataValue(NDV)
        DataSet.FlushCache()

        DataSet = None

def get_netcdf_info_global(filename,var_name):

    nc_file = gdal.Open(filename)

    if nc_file is None:
        print "Failed to open file, check directory and try again."
        sys.exit()

    #There are more than two variables, so specifying the lwe_thickness variable

    if nc_file.GetSubDatasets() > 1:
        subdataset = 'NETCDF:"'+filename+'":'+var_name #Specifying the subset name
        src_ds_sd = gdal.Open(subdataset) #Reading the subset
        NDV = src_ds_sd.GetRasterBand(1).GetNoDataValue() #Get the nodatavalues
        xsize = src_ds_sd.RasterXSize #Get the X size
        ysize = src_ds_sd.RasterYSize #Get the Y size
        GeoT = src_ds_sd.GetGeoTransform() #Get the GeoTransform
        Projection = osr.SpatialReference() #Get the SpatialReference
        Projection.ImportFromWkt(src_ds_sd.GetProjectionRef()) #Setting the Spatial Reference
        src_ds_sd = None #Closing the file
        nc_file = None #Closing the file

        return xsize,ysize,GeoT,NDV #Return data that will be used to convert the shapefile

#Upload GeoTiffs to geoserver
def upload_global_tiff(dir,geoserver_rest_url,workspace,uname,pwd):

    headers = {
        'Content-type': 'image/tiff',
    }

    dir = os.path.join(dir, '')
    for file in sorted(os.listdir(dir)): #Looping through all the files in the given directory
        if file is None:
            print "No files. Please check directory and try again."
            sys.exit()
        data = open(dir+file,'rb').read() #Read the file
        store_name = file.split('.')[0]  #Creating the store name dynamically
        request_url = '{0}workspaces/{1}/coveragestores/{2}/file.geotiff'.format(geoserver_rest_url,workspace,store_name) #Creating the rest url
        requests.put(request_url,verify=False,headers=headers,data=data,auth=(uname,pwd)) #Creating the resource on the geoserver

def finditem(obj, key):
    if key in obj: return obj[key]
    for k, v in obj.items():
        if isinstance(v,dict):
            return finditem(v, key)

def get_pt_plot(pt_coords):
    graph_json= {}

    ts_plot = []

    nc_file = '/grace/nepal/nepal.nc'

    coords = pt_coords.split(',')
    stn_lat = float(coords[1])
    stn_lon = float(coords[0])

    nc_fid = Dataset(nc_file,'r')
    nc_var = nc_fid.variables  # Get the netCDF variables
    nc_var.keys()  # Getting variable keys

    time = nc_var['time'][:]
    start_date = '01/01/2002'
    date_str = datetime.strptime(start_date, "%m/%d/%Y")  # Start Date string.
    lat = nc_var['lat'][:]
    lon = nc_var['lon'][:]

    for timestep, v in enumerate(time):

        current_time_step = nc_var['lwe_thickness'][timestep, :, :]  # Getting the index of the current timestep

        end_date = date_str + timedelta(days=float(v))  # Actual human readable date of the timestep

        data = nc_var['lwe_thickness'][timestep,:,:]

        lon_idx = (np.abs(lon - stn_lon)).argmin()
        lat_idx = (np.abs(lat - stn_lat)).argmin()

        value = data[lat_idx,lon_idx]

        time_stamp = calendar.timegm(end_date.utctimetuple()) * 1000

        ts_plot.append([time_stamp,round(float(value),3)])
        ts_plot.sort()

    graph_json["values"] = ts_plot
    graph_json["point"] = [round(stn_lat,2),round(stn_lon,2)]
    graph_json = json.dumps(graph_json)

    return graph_json


def get_global_plot(pt_coords):
    graph_json = {}

    ts_plot = []

    nc_file = get_global_nc()

    coords = pt_coords.split(',')
    stn_lat = float(coords[1])
    stn_lon = float(coords[0])

    nc_fid = Dataset(nc_file, 'r')
    nc_var = nc_fid.variables  # Get the netCDF variables
    nc_var.keys()  # Getting variable keys


    time = nc_var['time'][:]
    start_date = '01/01/2002'
    date_str = datetime.strptime(start_date, "%m/%d/%Y")  # Start Date string.
    lat = nc_var['lat'][:]
    lon = nc_var['lon'][:]

    for timestep, v in enumerate(time):
        current_time_step = nc_var['lwe_thickness'][timestep, :, :]  # Getting the index of the current timestep

        end_date = date_str + timedelta(days=float(v))  # Actual human readable date of the timestep

        data = nc_var['lwe_thickness'][timestep, :, :]

        lon_idx = (np.abs(lon - stn_lon)).argmin()
        lat_idx = (np.abs(lat - stn_lat)).argmin()

        value = data[lat_idx, lon_idx]

        time_stamp = calendar.timegm(end_date.utctimetuple()) * 1000

        ts_plot.append([time_stamp, round(float(value), 3)])
        ts_plot.sort()

    graph_json["values"] = ts_plot
    graph_json["point"] = [round(stn_lat, 2), round(stn_lon, 2)]
    graph_json = json.dumps(graph_json)
    return graph_json

def get_global_poly(bounds):
    graph_json = {}

    ts_plot = []

    nc_file = get_global_nc()

    miny = float(bounds[1])
    minx = float(bounds[0])
    maxx = float(bounds[2])
    maxy = float(bounds[3])

    nc_fid = Dataset(nc_file, 'r')
    nc_var = nc_fid.variables  # Get the netCDF variables
    nc_var.keys()  # Getting variable keys


    time = nc_var['time'][:]
    start_date = '01/01/2002'
    date_str = datetime.strptime(start_date, "%m/%d/%Y")  # Start Date string.
    lat = nc_var['lat'][:]
    lon = nc_var['lon'][:]

    for timestep, v in enumerate(time):
        current_time_step = nc_var['lwe_thickness'][timestep, :, :]  # Getting the index of the current timestep

        end_date = date_str + timedelta(days=float(v))  # Actual human readable date of the timestep

        data = nc_var['lwe_thickness'][timestep, :, :]

        lon_idx = (np.abs(lon - miny)).argmin()
        lat_idx = (np.abs(lat - minx)).argmin()
        lon2_idx = (np.abs(lon - maxy)).argmin()
        lat2_idx = (np.abs(lat - maxx)).argmin()

        values = data[lat_idx:lat2_idx,lon_idx:lon2_idx]

        if values is not None:

            value = np.mean(values)

            time_stamp = calendar.timegm(end_date.utctimetuple()) * 1000
            ts_plot.append([time_stamp, round(float(value), 3)])
            ts_plot.sort()

    ts_plot.sort()
    graph_json["values"] = ts_plot
    graph_json["point"] = [round(minx, 2), round(miny, 2),round(maxx, 2), round(maxy, 2)]
    graph_json = json.dumps(graph_json)
    return graph_json

def get_global_plot_api(pt_coords,start_date,end_date):
    graph_json = {}

    ts_plot = []

    nc_file = get_global_nc()

    coords = pt_coords.split(',')
    stn_lat = float(coords[1])
    stn_lon = float(coords[0])

    nc_fid = Dataset(nc_file, 'r')
    nc_var = nc_fid.variables  # Get the netCDF variables
    nc_var.keys()  # Getting variable keys

    time = nc_var['time'][:]
    start_date = '2002-01-01'
    date_str = datetime.strptime(start_date, "%Y-%m-%d")  # Start Date string.
    lat = nc_var['lat'][:]
    lon = nc_var['lon'][:]

    for timestep, v in enumerate(time):
        current_time_step = nc_var['lwe_thickness'][timestep, :, :]  # Getting the index of the current timestep

        actual_date = date_str + timedelta(days=float(v))  # Actual human readable date of the timestep


        data = nc_var['lwe_thickness'][timestep, :, :]

        lon_idx = (np.abs(lon - stn_lon)).argmin()
        lat_idx = (np.abs(lat - stn_lat)).argmin()

        value = data[lat_idx, lon_idx]

        time_stamp = calendar.timegm(actual_date.utctimetuple()) * 1000
        if start_date < unicode(actual_date) < end_date:
            ts_plot.append([time_stamp, round(float(value), 3)])
            ts_plot.sort()

    graph_json["values"] = ts_plot
    graph_json["point"] = [round(stn_lat, 2), round(stn_lon, 2)]
    graph_json = json.dumps(graph_json)
    return graph_json

def get_color_bar():

    value_range = [-50,50]
    min = value_range[0]
    max = value_range[1]
    opacity = [0.7] * 21

    cbar = ["#67001f",
            "#850c1e",
            "#a3201d",
            "#bd361c",
            "#d2501d",
            "#df6e22",
            "#e88e30",
            "#f0aa49",
            "#f7c670",
            "#fde1a6",
            "#fafafa",
            "#b7edf8",
            "#91d8f8",
            "#74bff9",
            "#5ea6f9",
            "#498dfa",
            "#3172fa",
            "#175be9",
            "#114ac0",
            "#0c3c94",
            "#053061"]

    interval = abs(min/10)

    scale = [x for x in range(min,max+1,interval)]
    final_cbar = zip(cbar,scale,opacity)

    return final_cbar

def get_pt_region(pt_coords,nc_file):

    graph_json = {}
    ts_plot = []

    coords = pt_coords.split(',')
    stn_lat = float(coords[1])
    stn_lon = float(coords[0])

    nc_fid = Dataset(nc_file, 'r')
    nc_var = nc_fid.variables  # Get the netCDF variables
    nc_var.keys()  # Getting variable keys

    time = nc_var['time'][:]
    start_date = '01/01/2002'
    date_str = datetime.strptime(start_date, "%m/%d/%Y")  # Start Date string.
    lat = nc_var['lat'][:]
    lon = nc_var['lon'][:]

    for timestep, v in enumerate(time):
        current_time_step = nc_var['lwe_thickness'][timestep, :, :]  # Getting the index of the current timestep

        end_date = date_str + timedelta(days=float(v))  # Actual human readable date of the timestep

        data = nc_var['lwe_thickness'][timestep, :, :]

        lon_idx = (np.abs(lon - stn_lon)).argmin()
        lat_idx = (np.abs(lat - stn_lat)).argmin()

        value = data[lat_idx, lon_idx]

        time_stamp = calendar.timegm(end_date.utctimetuple()) * 1000

        ts_plot.append([time_stamp, round(float(value), 3)])
        ts_plot.sort()

    graph_json["values"] = ts_plot
    graph_json["point"] = [round(stn_lat, 2), round(stn_lon, 2)]
    graph_json = json.dumps(graph_json)


    return graph_json

#Convert the shapefiles into a geojson object
def convert_shp(files):

    #Initizalizing an empty geojson string.
    geojson_string = ''

    try:
        #Storing the uploaded files in a temporary directory
        temp_dir = tempfile.mkdtemp()
        for f in files:
            f_name = f.name
            f_path = os.path.join(temp_dir,f_name)

            with open(f_path,'wb') as f_local:
                f_local.write(f.read())

        #Getting a list of files within the temporary directory
        for file in os.listdir(temp_dir):
            #Reading the shapefile only
            if file.endswith(".shp"):
                f_path = os.path.join(temp_dir,file)
                omit = ['SHAPE_AREA', 'SHAPE_LEN']

                #Reading the shapefile with fiona and reprojecting it
                with fiona.open(f_path) as source:
                    project = functools.partial(pyproj.transform,
                                                pyproj.Proj(**source.crs),
                                                pyproj.Proj(init='epsg:3857'))
                    features = []
                    for f in source:
                        shape = shapely.geometry.shape(f['geometry']) #Getting the shape of the shapefile
                        projected_shape = shapely.ops.transform(project, shape) #Transforming the shapefile

                        # Remove the properties we don't want
                        props = f['properties']  # props is a reference
                        for k in omit:
                            if k in props:
                                del props[k]

                        feature = geojson.Feature(id=f['id'],
                                                  geometry=projected_shape,
                                                  properties=props) #Creating a geojson feature by extracting properties through the fiona and shapely.geometry module
                        features.append(feature)
                    fc = geojson.FeatureCollection(features)

                    geojson_string = geojson.dumps(fc) #Creating the geojson string


    except:
        return 'error'
    finally:
        #Delete the temporary directory once the geojson string is created
        if temp_dir is not None:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

    return geojson_string

def convert_shp_bounds(bounds):
    inProj = Proj(init='epsg:3857')
    outProj = Proj(init='epsg:4326')
    minx, miny = bounds[0], bounds[1]
    maxx, maxy = bounds[2], bounds[3]
    x1, y1 = transform(inProj, outProj, minx, miny)
    x2, y2 = transform(inProj, outProj, maxx, maxy)

    x1 = x1 + 360
    x2 = x2 + 360
    reproj_bounds = [x1,y1,x2,y2]

    return reproj_bounds

def get_global_dates():
    grace_layer_options = []
    grace_nc = None
    for file in os.listdir(GLOBAL_NETCDF_DIR):
        if file.startswith('GRC') and file.endswith('.nc'):
            grace_nc = GLOBAL_NETCDF_DIR + file

    start_date = '01/01/2002'

    nc_fid = Dataset(grace_nc, 'r')  # Reading the netcdf file
    nc_var = nc_fid.variables  # Get the netCDF variables
    nc_var.keys()  # Getting variable keys

    time = nc_var['time'][:]

    date_str = datetime.strptime(start_date, "%m/%d/%Y")  # Start Date string.

    for timestep, v in enumerate(time):
        current_time_step = nc_var['lwe_thickness'][timestep, :, :]  # Getting the index of the current timestep

        end_date = date_str + timedelta(days=float(v))  # Actual human readable date of the timestep

        ts_file_name = end_date.strftime("%Y_%m_%d")  # Changing the date string format
        ts_display = end_date.strftime("%Y %B %d")
        grace_layer_options.append([ts_display,ts_file_name])

    return grace_layer_options
