import rasterio
import fiona
# from rasterio.tools.mask import mask
from netCDF4 import Dataset
import netCDF4, math,rtree
import os
from datetime import datetime,timedelta
import time
import calendar
import numpy as np
import shapefile as sf
import tempfile, shutil,sys
import gdal
import ogr
import osr
import requests
import csv, json
from grace import get_netcdf_info
from .app import Grace
import json
import functools
import shapely.geometry #Need this to find the bounds of a given geometry
import shapely.ops
import geojson
import pyproj
from pyproj import Proj, transform
from matplotlib.path import Path
from django.http import JsonResponse, HttpResponse, Http404

#Check if the user is superuser or staff. Only the superuser or staff have the permission to add and manage watersheds.
def user_permission_test(user):
    return user.is_superuser or user.is_staff
#
# def create_global_geotiffs(file_dir,geotiff_dir):
#
#     # Specify the relative file location
#     start_date = '01/01/2002'  # Date that GRACE data is available from
#
#     for file in os.listdir(file_dir): #Looping through the directory
#         if file is None:
#             print "No files to parse"
#             sys.exit()
#         nc_fid = Dataset(file_dir+file,'r') #Reading the netcdf file
#         nc_var = nc_fid.variables #Get the netCDF variables
#         nc_var.keys() #Getting variable keys
#         time = nc_var['time'][:] #Get the all the avaialable timesteps. Timestep increment value is x days after startdate
#
#         lwe_thickness = nc_var['lwe_thickness'][:,:,:] #Array with the all the values for lwe_thickness
#
#         date_str = datetime.datetime.strptime(start_date, "%m/%d/%Y") #Start Date string.
#
#
#         var = "lwe_thickness" #Specifying the variable key. This parameter will be used to retrieve information about the netCDF file
#         xsize, ysize, GeoT, Projection, NDV = get_netcdf_info(file_dir+file, var) #Get information about the netCDF file
#
#         ts_one = nc_var['lwe_thickness'][0,:,:]
#
#
#
#
#         unique_vals = set()
#         for i in ts_one:
#             for j in i:
#                 if float(j) not in unique_vals:
#                     unique_vals.add(float(j))
#
#         x = []
#         y = []
#         for i in unique_vals:
#             idx = np.where(nc_var['lwe_thickness'][0,:,:] == float(i))
#             x = x + idx[0].tolist()
#             y = y + idx[1].tolist()
#
#         x_y = zip(x,y)
#
#         grace_points = []
#         for i in x_y:
#             grace_json = {}  # Empty json object to store the corresponding latitude, longitude and lwe thickness value
#             latitude = nc_var['lat'][i[0]]
#             longitude = nc_var['lon'][i[1]]
#             thickness = nc_var['lwe_thickness'][0, i[0], i[1]]
#
#             # Saving all the values to the jon dictionary
#             grace_json["latitude"] = latitude
#             grace_json["longitude"] = longitude
#             grace_json["thickness"] = thickness
#             grace_points.append(grace_json)
#
#             # Creating the shapefile from the json dictionaries, then converting it to a raster
#         try:
#             file_name = 'grace_sites'
#             temp_dir = tempfile.mkdtemp()  # Creating a temporary directory to save the shapefile
#             file_location = temp_dir + "/" + file_name
#
#             w = sf.Writer(sf.POINT)  # Creating a point shapefile
#             w.field('thickness')  # Creating an attribute field called thickness for storing the variable value
#
#             # Looping through the list of json dictionaries to create points
#             for item in grace_points:
#                 w.point(float(item['longitude']), float(item['latitude']))  # Creating the point
#                 w.record(item['thickness'], 'Point')  # Assigning value to the point
#             w.save(file_location)
#
#             # Creating a projection file for the shapefile
#             prj = open("%s.prj" % file_location, "w")
#             epsg = 'GEOGCS["WGS84",DATUM["WGS_1984",SPHEROID["WGS84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'
#             prj.write(epsg)
#             prj.close()
#
#             # Begin the conveersion to a raster
#
#             NoData_value = -9999  # Specifying no data value
#             shp_file = file_location + ".shp"  # Find the shapefile location
#
#             out_loc = geotiff_dir + "tsone" + ".tif"  # Specify the GeoTiff name and output
#
#             source_ds = ogr.Open(shp_file)  # Reading the shapefile
#             source_layer = source_ds.GetLayer()  # Getting the actual layer
#             spatialRef = source_layer.GetSpatialRef()  # Get the Spatial Reference
#
#             raster_layer = gdal.GetDriverByName('GTiff').Create(out_loc, xsize, ysize, 1,
#                                                                 gdal.GDT_Float32)  # Initializing an empty GeoTiff
#             raster_layer.SetProjection(
#                 spatialRef.ExportToWkt())  # Set the projection based on the shapefile projection
#             raster_layer.SetGeoTransform(GeoT)  # Set the Geotransform.
#             band = raster_layer.GetRasterBand(1)  # Specifying the number of bands
#             band.SetNoDataValue(NoData_value)  # Setting no data values
#
#             band.FlushCache()  # This call will recover memory used to cache data blocks for this raster band, and ensure that new requests are referred to the underlying driver.
#
#             gdal.RasterizeLayer(raster_layer, [1], source_layer,
#                                 options=["ATTRIBUTE=thickness"])  # Create the GeoTiff layer
#         except:
#             print "Error parsing the data. Please check directory and try again."
#             sys.exit()
#             return False

def create_global_tiff(file_name,output_dir,var_name):

    # NewFile = '/home/tethys/geotiff_global/'+var_name+'.tif'
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
def upload_global_tiff(dir,geoserver_rest_url,workspace):

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
        requests.put(request_url,headers=headers,data=data,auth=('admin','geoserver')) #Creating the resource on the geoserver

# def clip_raster():
#     geoms = [{
#         "type": "Polygon",
#         "coordinates": [
#             [
#                 [
#                     79.78271484375,
#                     25.997549919572112
#                 ],
#                 [
#                     88.57177734375,
#                     25.997549919572112
#                 ],
#                 [
#                     88.57177734375,
#                     30.713503990354965
#                 ],
#                 [
#                     79.78271484375,
#                     30.713503990354965
#                 ],
#                 [
#                     79.78271484375,
#                     25.997549919572112
#                 ]
#             ]
#         ]
#     }]
#
#     file_input_dir = '/home/tethys/geotiff/'
#     file_outut_dir = '/home/tethys/geotiff_clipped/'
#     for file in os.listdir(file_input_dir):
#         with rasterio.open(file_input_dir+file) as src:
#             out_image, out_transform = mask(src, geoms, crop=True)
#         out_meta = src.meta.copy()
#         out_meta.update({"driver": "GTiff",
#                          "height": out_image.shape[1],
#                          "width": out_image.shape[2],
#                          "transform": out_transform})
#
#         with rasterio.open(file_outut_dir+file, "w", **out_meta) as dest:
#             dest.write(out_image)
#
# def clip_world():
#     file_input_dir = '/home/tethys/geotiff_global/'
#     file_outut_dir = '/home/tethys/geotiff_clipped_global/'
#
#     with fiona.open("/home/tethys/Downloads/world/TM_WORLD_BORDERS_SIMPL-0.3.shp","r") as shpfile:
#         features = [feature["geometry"] for feature in shpfile]
#
#     for file in os.listdir(file_input_dir):
#         with rasterio.open(file_input_dir+file) as src:
#             out_image, out_transform = mask(src, features, crop=True)
#         out_meta = src.meta.copy()
#         out_meta.update({"driver": "GTiff",
#                          "height": out_image.shape[1],
#                          "width": out_image.shape[2],
#                          "transform": out_transform})
#
#         with rasterio.open(file_outut_dir+file, "w", **out_meta) as dest:
#             dest.write(out_image)

def create_world_json():
    # file_input_dir = '/home/tethys/geotiff_global/'
    # file_outut_dir = '/home/tethys/geotiff_clipped_global/'

    with open('/home/tethys/Downloads/world.geojson') as f:
        data = json.load(f)

    for feature in data['features']:
        coordinates = feature['geometry']['coordinates']
        type = feature['geometry']['type']
        if type == 'Polygon':
            # coordinate = coordinates[0][0]
            for elem in coordinates:
                for coord in elem:
                    if coord[0] < 0.25:
                        coord[0] = coord[0]+360
        if type == 'MultiPolygon':
            for polygon in coordinates:
                for elem in polygon:
                    for coord in elem:
                        if coord[0] < 0.25:
                            coord[0] = coord[0]+360

                            # features = [feature['geometry'] for feature in data['features']]
                            # for file in os.listdir(file_input_dir):
                            #     with rasterio.open(file_input_dir+file) as src:
                            #         out_image, out_transform = mask(src, features, crop=True)
                            #     out_meta = src.meta.copy()
                            #     out_meta.update({"driver": "GTiff",
                            #                      "height": out_image.shape[1],
                            #                      "width": out_image.shape[2],
                            #                      "transform": out_transform})
                            #     with rasterio.open(file_outut_dir+file, "w", **out_meta) as dest:
                            #         dest.write(out_image)

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

    nc_file = '/grace/global/GRCTellus.JPL.200204_201608.GLO.RL05M_1.MSCNv02CRIv02.nc'

    # nc_file = '/home/tethys/netcdf/global/GRCTellus.JPL.200204_201608.GLO.RL05M_1.MSCNv02CRIv02.nc'

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

    nc_file = '/grace/global/GRCTellus.JPL.200204_201608.GLO.RL05M_1.MSCNv02CRIv02.nc'

    # nc_file = '/home/tethys/netcdf/global/GRCTellus.JPL.200204_201608.GLO.RL05M_1.MSCNv02CRIv02.nc'

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

    nc_file = '/grace/global/GRCTellus.JPL.200204_201608.GLO.RL05M_1.MSCNv02CRIv02.nc'

    # nc_file = '/home/tethys/netcdf/global/GRCTellus.JPL.200204_201608.GLO.RL05M_1.MSCNv02CRIv02.nc'
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


def process_shp(gbyos_grc_ncf,gbyos_fct_ncf,bounds):
    graph_json = {}
    ts_plot = []

    f = netCDF4.Dataset(gbyos_grc_ncf, 'r')

    IS_grc_lon = len(f.dimensions['lon'])

    IS_grc_lat = len(f.dimensions['lat'])

    IS_grc_time = len(f.dimensions['time'])

    ZV_grc_lon = f.variables['lon']
    ZV_grc_lat = f.variables['lat']
    ZV_grc_time = f.variables['time']

    ZS_grc_fil = netCDF4.default_fillvals['f4']
    if 'RUNSF' in f.variables:
        var = f.variables['RUNSF']
        if '_FillValue' in var.ncattrs():
            ZS_grc_fil = var._FillValue
            print(' - The fill value for RUNSF is: ' + str(ZS_grc_fil))
        else:
            ZS_grc_fil = None

    g = netCDF4.Dataset(gbyos_fct_ncf, 'r')
    IS_fct_lon = len(g.dimensions['lon'])
    IS_fct_lat = len(g.dimensions['lat'])
    ZV_fct_lon = g.variables['lon']
    ZV_fct_lat = g.variables['lat']

    index = rtree.index.Index()
    shp_bounds = []
    miny = float(bounds[1])
    minx = float(bounds[0])
    maxx = float(bounds[2])
    maxy = float(bounds[3])
    bbPath = Path(np.array([[minx, miny],
                            [maxx, miny],
                            [maxx, maxy],
                            [minx, maxy],
                            ]))

    gbyos_pol_fid = int(1)
    index.insert(gbyos_pol_fid, [minx,miny,maxx,maxy])
    shp_bounds.append(bounds)
    IS_dom_tot = 0
    IV_dom_lon = []
    IV_dom_lat = []

    for JS_grc_lon in range(IS_grc_lon):
        ZS_grc_lon = ZV_grc_lon[JS_grc_lon]
        if (ZS_grc_lon > 180):
            ZS_grc_lon = ZS_grc_lon - 360
            # Shift GRACE longitude range from [0;360] to [-180;180]

        for JS_grc_lat in range(IS_grc_lat):
            ZS_grc_lat = ZV_grc_lat[JS_grc_lat]
            gbyos_pnt_shy = shapely.geometry.Point(ZS_grc_lon, ZS_grc_lat)
            # a shapely point now exists for a given GRACE grid cell
            for gbyos_pol_fid in [int(x) for x in list(index.intersection(gbyos_pnt_shy.bounds))]:
                gbyos_pol_fea = gbyos_pol_fid
                if (bbPath.contains_point((float(ZS_grc_lon), float(ZS_grc_lat)))):
                    IV_dom_lon.append(JS_grc_lon)
                    IV_dom_lat.append(JS_grc_lat)
                    IS_dom_tot = IS_dom_tot + 1

    print(' - The number of grid cells found is: ' + str(IS_dom_tot))
    ZV_dom_avg = [0] * IS_dom_tot
    for JS_dom_tot in range(IS_dom_tot):
        JS_grc_lon = IV_dom_lon[JS_dom_tot]
        JS_grc_lat = IV_dom_lat[JS_dom_tot]
        for JS_grc_time in range(IS_grc_time):
            ZV_dom_avg[JS_dom_tot] = ZV_dom_avg[JS_dom_tot] \
                                     + f.variables['lwe_thickness'] \
                                         [JS_grc_time, JS_grc_lat, JS_grc_lon]
    ZV_dom_avg = [x / IS_grc_time for x in ZV_dom_avg]

    ZV_dom_sqm = [0] * IS_dom_tot
    for JS_dom_tot in range(IS_dom_tot):
        JS_grc_lat = IV_dom_lat[JS_dom_tot]
        ZS_grc_lat = ZV_grc_lat[JS_grc_lat]
        ZV_dom_sqm[JS_dom_tot] = 6371000 * math.radians(0.5) \
                                 * 6371000 * math.radians(0.5) \
                                 * math.cos(math.radians(ZS_grc_lat))
    ZM_grc_scl = g.variables['scale_factor'][:, :]
    IS_dom_msk = 0
    ZS_sqm = 0
    for JS_dom_tot in range(IS_dom_tot):
        JS_grc_lon = IV_dom_lon[JS_dom_tot]
        JS_grc_lat = IV_dom_lat[JS_dom_tot]
        if (ZM_grc_scl.mask[JS_grc_lat, JS_grc_lon]):
            IS_dom_msk = IS_dom_msk + 1
        else:
            ZS_sqm = ZS_sqm + ZV_dom_sqm[JS_dom_tot]

    ZV_wsa = []
    for JS_grc_time in range(IS_grc_time):
        ZS_wsa = 0
        for JS_dom_tot in range(IS_dom_tot):
            JS_grc_lon = IV_dom_lon[JS_dom_tot]
            JS_grc_lat = IV_dom_lat[JS_dom_tot]
            ZS_dom_sqm = ZV_dom_sqm[JS_dom_tot]
            ZS_dom_avg = ZV_dom_avg[JS_dom_tot]
            if (ZM_grc_scl.mask[JS_grc_lat, JS_grc_lon]):
                ZS_dom_scl = 0
            else:
                ZS_dom_scl = ZM_grc_scl[JS_grc_lat, JS_grc_lon]
            ZS_dom_wsa = (f.variables['lwe_thickness'] \
                              [JS_grc_time, JS_grc_lat, JS_grc_lon] \
                          - ZS_dom_avg) / 100 \
                         * ZS_dom_scl \
                         * ZS_dom_sqm
            # The division by 100 is to go from cm to m in GRACE data.
            ZS_wsa = ZS_wsa + ZS_dom_wsa
        ZV_wsa.append(100 * ZS_wsa / ZS_sqm)

    gbyos_dat_str = datetime.strptime('2002-01-01T00:00:00', '%Y-%m-%dT%H:%M:%S')

    YV_grc_time = []
    for JS_grc_time in range(IS_grc_time):
        gbyos_dat_dlt = timedelta(days=ZV_grc_time[JS_grc_time])
        YS_grc_time = (gbyos_dat_str + gbyos_dat_dlt).strftime('%m/%d/%Y')
        YV_grc_time.append(YS_grc_time)

    for JS_grc_time in range(IS_grc_time):
        time_stamp = time.mktime(datetime.strptime(YV_grc_time[JS_grc_time], "%m/%d/%Y").timetuple()) * 1000
        ts_plot.append([time_stamp, round(float(ZV_wsa[JS_grc_time]), 3)])



    ts_plot.sort()
    graph_json["values"] = ts_plot
    graph_json["point"] = [round(minx, 2), round(miny, 2), round(maxx, 2), round(maxy, 2)]
    graph_json = json.dumps(graph_json)
    return graph_json

def vals_from_shp(shapefile,GLOBAL_NETCDF_DIR):

    try:
        GLOBAL_NETCDF_DIR = os.path.join(GLOBAL_NETCDF_DIR, '')

        temp_dir = tempfile.mkdtemp()
        for f in shapefile:
            f_name = f.name
            f_path = os.path.join(temp_dir, f_name)

            with open(f_path, 'wb') as f_local:
                f_local.write(f.read())

        for file in os.listdir(temp_dir):
            # Reading the shapefile only
            if file.endswith(".shp"):
                f_path = os.path.join(temp_dir, file)
                gbyos_pol_shp = f_path

        for file in os.listdir(GLOBAL_NETCDF_DIR):
            if file.startswith('GRC') and file.endswith('.nc'):
                gbyos_grc_ncf = GLOBAL_NETCDF_DIR + file
            if file.startswith('CLM4') and file.endswith('.nc'):
                gbyos_fct_ncf = GLOBAL_NETCDF_DIR + file


        # *******************************************************************************
        # Read GRACE netCDF file
        # *******************************************************************************

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Open netCDF file
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        print('Read GRACE netCDF file')
        f = netCDF4.Dataset(gbyos_grc_ncf, 'r')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Get dimension sizes
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IS_grc_lon = len(f.dimensions['lon'])
        print(' - The number of longitudes is: ' + str(IS_grc_lon))
        IS_grc_lat = len(f.dimensions['lat'])
        print(' - The number of latitudes is: ' + str(IS_grc_lat))
        IS_grc_time = len(f.dimensions['time'])
        print(' - The number of time steps is: ' + str(IS_grc_time))
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Get values of dimension arrays
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ZV_grc_lon = f.variables['lon']
        ZV_grc_lat = f.variables['lat']
        ZV_grc_time = f.variables['time']

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Get fill values
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ZS_grc_fil = netCDF4.default_fillvals['f4']
        if 'RUNSF' in f.variables:
            var = f.variables['RUNSF']
            if '_FillValue' in var.ncattrs():
                ZS_grc_fil = var._FillValue
                print(' - The fill value for RUNSF is: ' + str(ZS_grc_fil))
            else:
                ZS_grc_fil = None

        # *******************************************************************************
        # Read scale factors netCDF file
        # *******************************************************************************

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Open netCDF file
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        print('Read scale factors netCDF file')
        g = netCDF4.Dataset(gbyos_fct_ncf, 'r')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Get dimension sizes
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IS_fct_lon = len(g.dimensions['lon'])
        print(' - The number of longitudes is: ' + str(IS_fct_lon))
        IS_fct_lat = len(g.dimensions['lat'])
        print(' - The number of latitudes is: ' + str(IS_fct_lat))
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Get values of dimension arrays
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ZV_fct_lon = g.variables['lon']
        ZV_fct_lat = g.variables['lat']

        # *******************************************************************************
        # Read polygon shapefile
        # *******************************************************************************
        print('Read polygon shapefile')
        gbyos_pol_lay = fiona.open(gbyos_pol_shp, 'r')
        IS_pol_tot = len(gbyos_pol_lay)
        print(' - The number of polygon features is: ' + str(IS_pol_tot))
        # *******************************************************************************
        # Create spatial index for the bounds of each polygon feature
        # *******************************************************************************


        index = rtree.index.Index()
        shp_bounds = []

        def explode(coords):
            """Explode a GeoJSON geometry's coordinates object and yield coordinate tuples.
            As long as the input is conforming, the type of the geometry doesn't matter."""
            for e in coords:
                if isinstance(e, (float, int, long)):
                    yield coords
                    break
                else:
                    for f in explode(e):
                        yield f

        def bbox(f):
            x, y = zip(*list(explode(f['geometry']['coordinates'])))
            return min(x), min(y), max(x), max(y)

        for gbyos_pol_fea in gbyos_pol_lay:
            gbyos_pol_fid = int(gbyos_pol_fea['id'])
            # the first argument of index.insert has to be 'int', not 'long' or 'str'
            gbyos_pol_shy = shapely.geometry.shape(gbyos_pol_fea['geometry'])
            index.insert(gbyos_pol_fid, gbyos_pol_shy.bounds)
            shp_bounds.append(gbyos_pol_shy.bounds)
            bbox_val = bbox(gbyos_pol_fea)
            # creates an index between the feature ID and the bounds of that feature

        # *******************************************************************************
        # Find GRACE grid cells that intersect with polygon
        # *******************************************************************************
        print('Find GRACE grid cells that intersect with polygon')
        IS_dom_tot = 0
        IV_dom_lon = []
        IV_dom_lat = []

        for JS_grc_lon in range(IS_grc_lon):
            ZS_grc_lon = ZV_grc_lon[JS_grc_lon]
            if (ZS_grc_lon > 180):
                ZS_grc_lon = ZS_grc_lon - 360
                # Shift GRACE longitude range from [0;360] to [-180;180]

            for JS_grc_lat in range(IS_grc_lat):
                ZS_grc_lat = ZV_grc_lat[JS_grc_lat]
                gbyos_pnt_shy = shapely.geometry.Point(ZS_grc_lon, ZS_grc_lat)
                # a shapely point now exists for a given GRACE grid cell
                for gbyos_pol_fid in \
                        [int(x) for x in list(index.intersection(gbyos_pnt_shy.bounds))]:
                    gbyos_pol_fea = gbyos_pol_lay[gbyos_pol_fid]
                    gbyos_pol_shy = shapely.geometry.shape(gbyos_pol_fea['geometry'])
                    if (gbyos_pnt_shy.within(gbyos_pol_shy)):
                        IV_dom_lon.append(JS_grc_lon)
                        IV_dom_lat.append(JS_grc_lat)
                        IS_dom_tot = IS_dom_tot + 1

        print(' - The number of grid cells found is: ' + str(IS_dom_tot))
        # *******************************************************************************
        # Find long-term mean for each intersecting GRACE grid cell
        # *******************************************************************************
        print('Find long-term mean for each intersecting GRACE grid cell')
        ZV_dom_avg = [0] * IS_dom_tot
        for JS_dom_tot in range(IS_dom_tot):
            JS_grc_lon = IV_dom_lon[JS_dom_tot]
            JS_grc_lat = IV_dom_lat[JS_dom_tot]
            for JS_grc_time in range(IS_grc_time):
                ZV_dom_avg[JS_dom_tot] = ZV_dom_avg[JS_dom_tot] \
                                         + f.variables['lwe_thickness'] \
                                             [JS_grc_time, JS_grc_lat, JS_grc_lon]
        ZV_dom_avg = [x / IS_grc_time for x in ZV_dom_avg]

        # *******************************************************************************
        # Compute surface area of each grid cell
        # *******************************************************************************
        print('Compute surface area of each grid cell')

        ZV_dom_sqm = [0] * IS_dom_tot
        for JS_dom_tot in range(IS_dom_tot):
            JS_grc_lat = IV_dom_lat[JS_dom_tot]
            ZS_grc_lat = ZV_grc_lat[JS_grc_lat]
            ZV_dom_sqm[JS_dom_tot] = 6371000 * math.radians(0.5) \
                                     * 6371000 * math.radians(0.5) \
                                     * math.cos(math.radians(ZS_grc_lat))

        # *******************************************************************************
        # Find number of NoData points in scale factors for shapefile and area
        # *******************************************************************************
        print('Find number of NoData points in scale factors for shapefile and area')
        ZM_grc_scl = g.variables['scale_factor'][:, :]
        IS_dom_msk = 0
        ZS_sqm = 0
        for JS_dom_tot in range(IS_dom_tot):
            JS_grc_lon = IV_dom_lon[JS_dom_tot]
            JS_grc_lat = IV_dom_lat[JS_dom_tot]
            if (ZM_grc_scl.mask[JS_grc_lat, JS_grc_lon]):
                IS_dom_msk = IS_dom_msk + 1
            else:
                ZS_sqm = ZS_sqm + ZV_dom_sqm[JS_dom_tot]
        print(' - The number of NoData points found is: ' + str(IS_dom_msk))
        print(' - The area (m2) for the domain is: ' + str(ZS_sqm))
        # *******************************************************************************
        # Compute total terrestrial water storage anomaly timeseries
        # *******************************************************************************
        print('Compute total terrestrial water storage anomaly timeseries')
        ZV_wsa = []
        for JS_grc_time in range(IS_grc_time):
            ZS_wsa = 0
            for JS_dom_tot in range(IS_dom_tot):
                JS_grc_lon = IV_dom_lon[JS_dom_tot]
                JS_grc_lat = IV_dom_lat[JS_dom_tot]
                ZS_dom_sqm = ZV_dom_sqm[JS_dom_tot]
                ZS_dom_avg = ZV_dom_avg[JS_dom_tot]
                if (ZM_grc_scl.mask[JS_grc_lat, JS_grc_lon]):
                    ZS_dom_scl = 0
                else:
                    ZS_dom_scl = ZM_grc_scl[JS_grc_lat, JS_grc_lon]
                ZS_dom_wsa = (f.variables['lwe_thickness'] \
                                  [JS_grc_time, JS_grc_lat, JS_grc_lon] \
                              - ZS_dom_avg) / 100 \
                             * ZS_dom_scl \
                             * ZS_dom_sqm
                # The division by 100 is to go from cm to m in GRACE data.
                ZS_wsa = ZS_wsa + ZS_dom_wsa
            ZV_wsa.append(100 * ZS_wsa / ZS_sqm)

        # *******************************************************************************
        # Determine time strings
        # *******************************************************************************
        print('Determine time strings')
        gbyos_dat_str = datetime.strptime('2002-01-01T00:00:00', '%Y-%m-%dT%H:%M:%S')

        YV_grc_time = []
        for JS_grc_time in range(IS_grc_time):
            gbyos_dat_dlt = timedelta(days=ZV_grc_time[JS_grc_time])
            YS_grc_time = (gbyos_dat_str + gbyos_dat_dlt).strftime('%m/%d/%Y')
            YV_grc_time.append(YS_grc_time)

        return gbyos_dat_str
    except Exception as e:

        if temp_dir is not None:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
        return JsonResponse({"error": e})

    finally:
        # Delete the temporary directory once the geojson string is created
        if temp_dir is not None:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

