from netCDF4 import Dataset
import os,os.path
from datetime import datetime,timedelta
import numpy as np
import shapefile as sf
import os, tempfile, shutil,sys
import gdal
import ogr
import osr
import requests, urlparse
import csv
from tethys_dataset_services.engines import GeoServerSpatialDatasetEngine
from model import *

def create_geotiffs(file_dir,geotiff_dir):

     #Specify the relative file location
    start_date = '01/01/2002' #Date that GRACE data is available from
    for file in os.listdir(file_dir): #Looping through the directory

        if file is None:
            print "No files to parse"
            sys.exit()
        if file.endswith('.nc'):
            nc_fid = Dataset(file_dir+file,'r') #Reading the netcdf file
            nc_var = nc_fid.variables #Get the netCDF variables
            nc_var.keys() #Getting variable keys
            time = nc_var['time'][:] #Get the all the avaialable timesteps. Timestep increment value is x days after startdate

            lwe_thickness = nc_var['lwe_thickness'][:,:,:] #Array with the all the values for lwe_thickness

            date_str = datetime.strptime(start_date, "%m/%d/%Y") #Start Date string.


            var = "lwe_thickness" #Specifying the variable key. This parameter will be used to retrieve information about the netCDF file
            xsize, ysize, GeoT, Projection, NDV = get_netcdf_info(file_dir+file, var) #Get information about the netCDF file

            #Looping through each timestep
            for timestep,v in enumerate(time):

                current_time_step = nc_var['lwe_thickness'][timestep, :, :] #Getting the index of the current timestep

                end_date = date_str + timedelta(days=float(v)) #Actual human readable date of the timestep

                ts_file_name = end_date.strftime("%Y_%m_%d") #Changing the date string format
                fts_vals = set() #Creating an empty set to store the values for the timestep

                for i in current_time_step.compressed(): #Compressed ignores null values and returns an array of values that exist
                    if float(i) not in fts_vals: #Check if the value exists in the fts_values. If not then add the values.
                        fts_vals.add(float(i))

                x = [] #Creating an empty list to store the x indexes for a given value
                y = [] #Creating an empty list to store the y indexes for a given value

                for i in fts_vals: #Looping through the existing values to get the indexes of x,y
                    idx = np.where(nc_var['lwe_thickness'][timestep, :, :] == float(i)) #Find the index of the given value
                    x = x + idx[0].tolist() #Write the x index values to the empty x list
                    y = y + idx[1].tolist() #Write the y index values to the empty y list

                x_y = zip(x, y) #Combining the x,y list

                grace_points = [] #Creating an empty list to store a list of json dictionaries. Will be used to generate the shapefile.


                for i in x_y: #Looping through the indexes to find the exact latitude and longitude of an existing value

                    grace_json = {} #Empty json object to store the corresponding latitude, longitude and lwe thickness value
                    latitude = nc_var['lat'][i[0]]
                    longitude = nc_var['lon'][i[1]]
                    thickness = nc_var['lwe_thickness'][timestep, i[0], i[1]]

                    #Saving all the values to the jon dictionary
                    grace_json["latitude"] = latitude
                    grace_json["longitude"] = longitude
                    grace_json["thickness"] = thickness
                    grace_points.append(grace_json)

                #Creating the shapefile from the json dictionaries, then converting it to a raster
                try:

                    file_name = 'grace_sites'
                    temp_dir = tempfile.mkdtemp() #Creating a temporary directory to save the shapefile
                    file_location = temp_dir + "/" + file_name

                    w = sf.Writer(sf.POINT) #Creating a point shapefile
                    w.field('thickness')  #Creating an attribute field called thickness for storing the variable value

                    #Looping through the list of json dictionaries to create points
                    for item in grace_points:
                        w.point(float(item['longitude']), float(item['latitude'])) #Creating the point
                        w.record(item['thickness'], 'Point') #Assigning value to the point
                    w.save(file_location)

                    #Creating a projection file for the shapefile
                    prj = open("%s.prj" % file_location, "w")
                    epsg = 'GEOGCS["WGS84",DATUM["WGS_1984",SPHEROID["WGS84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'
                    prj.write(epsg)
                    prj.close()

                    #Begin the conveersion to a raster

                    NoData_value = -9999 #Specifying no data value
                    shp_file = file_location + ".shp" #Find the shapefile location

                    out_loc = geotiff_dir+ts_file_name+ ".tif" #Specify the GeoTiff name and output

                    source_ds = ogr.Open(shp_file) #Reading the shapefile
                    source_layer = source_ds.GetLayer() #Getting the actual layer
                    spatialRef = source_layer.GetSpatialRef() #Get the Spatial Reference


                    raster_layer = gdal.GetDriverByName('GTiff').Create(out_loc, xsize, ysize, 1, gdal.GDT_Float32) #Initializing an empty GeoTiff
                    raster_layer.SetProjection(spatialRef.ExportToWkt()) #Set the projection based on the shapefile projection
                    raster_layer.SetGeoTransform(GeoT) #Set the Geotransform.
                    band = raster_layer.GetRasterBand(1) #Specifying the number of bands
                    band.SetNoDataValue(NoData_value) #Setting no data values

                    band.FlushCache() #This call will recover memory used to cache data blocks for this raster band, and ensure that new requests are referred to the underlying driver.

                    gdal.RasterizeLayer(raster_layer, [1], source_layer, options=["ATTRIBUTE=thickness"]) #Create the GeoTiff layer
                except:
                    print "Error parsing the data. Please check directory and try again."
                    sys.exit()
                    return False
                finally:
                    # Delete the temp shapefile dir after uploading the shapefile
                    if temp_dir is not None:
                        if os.path.exists(temp_dir):
                            shutil.rmtree(temp_dir)

#Get info from the netCDF file. This info will be used to convert the shapefile to a raster layer
def get_netcdf_info(filename,var_name):

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

        return xsize,ysize,GeoT,Projection,NDV #Return data that will be used to convert the shapefile

#Upload GeoTiffs to geoserver
def upload_tiff(dir,region,geoserver_rest_url,workspace,uname,pwd):
    print "just got to the upload tiff function"
    headers = {
        'Content-type': 'image/tiff',
    }
    spatial_data_engine = GeoServerSpatialDatasetEngine(endpoint=geoserver_rest_url,username=uname,password=pwd)

    #Check if workspace exists
    ws_name = workspace
    geoserver_uri = 'www.google.com'
    response = spatial_data_engine.list_workspaces()

    if response['success']:
        workspaces = response['result']

        if ws_name not in workspaces:
            spatial_data_engine.create_workspace(workspace_id=ws_name,uri=geoserver_uri)
    for file in os.listdir(dir): #Looping through all the files in the given directory
        if file is None:
            print "No files. Please check directory and try again."
            sys.exit()
        data = open(dir+file,'rb').read() #Read the file
        store_name = file.split('.')[0]+'_'+region #Creating the store name dynamically

        request_url = '{0}workspaces/{1}/coveragestores/{2}/file.geotiff'.format(geoserver_rest_url,ws_name,store_name) #Creating the rest url
        print request_url
        requests.put(request_url,headers=headers,data=data,auth=(uname,pwd)) #Creating the resource on the geoserver

def get_max_min(file_dir,output_dir):
    # Specify the relative file location
    start_date = '01/01/2002'  # Date that GRACE data is available from
    for file in os.listdir(file_dir):  # Looping through the directory
        if file is None:
            print "No files to parse"
            sys.exit()
        nc_fid = Dataset(file_dir + file, 'r')  # Reading the netcdf file
        nc_var = nc_fid.variables  # Get the netCDF variables
        nc_var.keys()  # Getting variable keys
        time = nc_var['time'][
               :]  # Get the all the avaialable timesteps. Timestep increment value is x days after startdate

        lwe_thickness = nc_var['lwe_thickness'][:, :, :]  # Array with the all the values for lwe_thickness

        date_str = datetime.strptime(start_date, "%m/%d/%Y")  # Start Date string.

        var = "lwe_thickness"  # Specifying the variable key. This parameter will be used to retrieve information about the netCDF file
        xsize, ysize, GeoT, Projection, NDV = get_netcdf_info(file_dir + file,
                                                              var)  # Get information about the netCDF file

        # Looping through each timestep
        with open(output_dir+"legend.csv","w") as f:
            writer = csv.writer(f)
            for timestep, v in enumerate(time):

                current_time_step = nc_var['lwe_thickness'][timestep, :, :]  # Getting the index of the current timestep

                end_date = date_str + timedelta(days=float(v))  # Actual human readable date of the timestep

                ts_file_name = end_date.strftime("%Y_%m_%d")  # Changing the date string format

                fts_vals = set()  # Creating an empty set to store the values for the timestep

                for i in current_time_step.compressed():  # Compressed ignores null values and returns an array of values that exist
                    if float(i) not in fts_vals:  # Check if the value exists in the fts_values. If not then add the values.
                        fts_vals.add(float(i))
                writer.writerow([ts_file_name,max(fts_vals),min(fts_vals)])