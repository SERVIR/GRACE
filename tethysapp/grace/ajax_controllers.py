from django.http import JsonResponse, HttpResponse, Http404
from django.shortcuts import render
from django.contrib.auth.decorators import login_required,user_passes_test
from django.views.decorators.csrf import csrf_exempt
from tethys_sdk.gizmos import *
from utilities import *
import json
from tethys_dataset_services.engines import GeoServerSpatialDatasetEngine
from sqlalchemy.orm.exc import ObjectDeletedError
from sqlalchemy.exc import IntegrityError
from model import *
import requests, urlparse
from gbyos import *
import shapely.geometry
import os

GRACE_NETCDF_DIR = '/grace/'
GLOBAL_NETCDF_DIR = '/grace/global/'

def plot_region(request):
    return_obj = {}
    if request.is_ajax() and request.method == 'POST':
        info = request.POST

        region_id = info.get('region-info')
        pt_coords = request.POST['point-lat-lon']

        Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
        session = Session()

        region = session.query(Region).get(region_id)
        display_name = region.display_name
        region_store = ''.join(display_name.split()).lower()

        FILE_DIR = os.path.join(GRACE_NETCDF_DIR, '')

        region_dir = os.path.join(FILE_DIR + region_store, '')

        nc_file = os.path.join(region_dir+region_store+".nc")

        if pt_coords:
            graph = get_pt_region(pt_coords,nc_file)
            graph = json.loads(graph)
            return_obj["values"] = graph["values"]
            return_obj["location"] = graph["point"]

        return_obj["success"] = "success"


    return JsonResponse(return_obj)


def get_plot(request):

    return_obj = {}

    if request.is_ajax() and request.method == 'POST':
        # Get the point/polygon/shapefile coordinates along with the selected variable
        pt_coords = request.POST['point-lat-lon']
        # poly_coords = request.POST['poly-lat-lon']
        # shp_bounds = request.POST['shp-lat-lon']

        if pt_coords:
            graph = get_pt_plot(pt_coords)
            graph = json.loads(graph)
            return_obj["values"] = graph["values"]
            return_obj["location"] = graph["point"]

        return_obj['success'] = "success"

    return JsonResponse(return_obj)

def get_plot_global(request):

    return_obj = {}

    if request.is_ajax() and request.method == 'POST':
        # Get the point/polygon/shapefile coordinates along with the selected variable
        pt_coords = request.POST['point-lat-lon']
        poly_coords = request.POST['poly-lat-lon']
        shp_bounds = request.POST['shp-lat-lon']

        GLOBAL_DIR = os.path.join(GLOBAL_NETCDF_DIR, '')


        for file in os.listdir(GLOBAL_DIR):
            if file.startswith('GRC') and file.endswith('.nc'):
                gbyos_grc_ncf = GLOBAL_NETCDF_DIR + file
            if file.startswith('CLM4') and file.endswith('.nc'):
                gbyos_fct_ncf = GLOBAL_NETCDF_DIR + file

        if pt_coords:
            graph = get_global_plot(pt_coords)
            graph = json.loads(graph)
            return_obj["values"] = graph["values"]
            return_obj["location"] = graph["point"]



        if poly_coords:
            poly_geojson = json.loads(poly_coords)
            # Get the bounds of the shape using shapely.geometry
            shape_obj = shapely.geometry.asShape(poly_geojson)
            poly_bounds = shape_obj.bounds
            graph = process_shp(gbyos_grc_ncf, gbyos_fct_ncf, poly_bounds)
            # graph = get_global_poly(poly_bounds)
            graph = json.loads(graph)
            return_obj["values"] = graph["values"]
            return_obj["location"] = graph["point"]

        if shp_bounds:
            shp_bounds = tuple(shp_bounds.split(','))
            # graph = get_global_poly(shp_bounds)

            graph = process_shp(gbyos_grc_ncf,gbyos_fct_ncf,shp_bounds)
            graph = json.loads(graph)
            return_obj["values"] = graph["values"]
            return_obj["location"] = graph["point"]

        return_obj['success'] = "success"

    return JsonResponse(return_obj)



@user_passes_test(user_permission_test)
def region_add(request):

    response = {}

    if request.is_ajax() and request.method == 'POST':
        info = request.POST

        region_name = info.get('region_name')
        region_store = ''.join(region_name.split()).lower()
        geoserver_id = info.get('geoserver')

        shapefile = request.FILES.getlist('shapefile')

        Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
        session = Session()


        geoserver = session.query(Geoserver).get(geoserver_id)
        url,uname,pwd = geoserver.url,geoserver.username,geoserver.password

        process_shapefile(shapefile, url, uname, pwd, region_store, GRACE_NETCDF_DIR, GLOBAL_NETCDF_DIR,region_name,geoserver_id)

        response = {"success":"success"}

        return JsonResponse(response)

@user_passes_test(user_permission_test)
def geoserver_add(request):

    response = {}

    if request.is_ajax() and request.method == 'POST':
        info = request.POST

        geoserver_name = info.get('geoserver_name')
        geoserver_url = info.get('geoserver_url')
        geoserver_username = info.get('geoserver_username')
        geoserver_password = info.get('geoserver_password')

        try:
            spatial_dataset_engine = GeoServerSpatialDatasetEngine(endpoint=geoserver_url,username=geoserver_username,password=geoserver_password)
            layer_list = spatial_dataset_engine.list_layers(debug=True)
            if layer_list:
                Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
                session = Session()
                geoserver = Geoserver(name=geoserver_name, url=geoserver_url, username=geoserver_username, password=geoserver_password)
                session.add(geoserver)
                session.commit()
                session.close()
                response = {"data": geoserver_name, "success": "Success"}
        except:
            response={"error":"Error processing the Geoserver URL. Please check the url,username and password."}


        return JsonResponse(response)


@user_passes_test(user_permission_test)
def geoserver_update(request):
    """
    Controller for updating a geoserver.
    """
    if request.is_ajax() and request.method == 'POST':
        # get/check information from AJAX request
        post_info = request.POST
        geoserver_id = post_info.get('geoserver_id')
        geoserver_name = post_info.get('geoserver_name')
        geoserver_url = post_info.get('geoserver_url')
        geoserver_username = post_info.get('geoserver_username')
        geoserver_password = post_info.get('geoserver_password')
        # check data
        if not geoserver_id or not geoserver_name or not geoserver_url or not \
                geoserver_username or not geoserver_password:
            return JsonResponse({'error': "Missing input data."})
        # make sure id is id
        try:
            int(geoserver_id)
        except ValueError:
            return JsonResponse({'error': 'Geoserver id is faulty.'})

        Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
        session = Session()

        geoserver = session.query(Geoserver).get(geoserver_id)
        try:
            spatial_dataset_engine = GeoServerSpatialDatasetEngine(endpoint=geoserver_url, username=geoserver_username,
                                                                   password=geoserver_password)
            layer_list = spatial_dataset_engine.list_layers(debug=True)
            if layer_list:


                geoserver.geoserver_name = geoserver_name
                geoserver.geoserver_url = geoserver_url
                geoserver.geoserver_username = geoserver_username
                geoserver.geoserver_password = geoserver_password

                session.commit()
                session.close()
                return JsonResponse({'success': "Geoserver sucessfully updated!"})
        except:
            return JsonResponse({'error': "A problem with your request exists."})


@user_passes_test(user_permission_test)
def geoserver_delete(request):
    """
    Controller for deleting a geoserver.
    """
    if request.is_ajax() and request.method == 'POST':
        # get/check information from AJAX request
        post_info = request.POST
        geoserver_id = post_info.get('geoserver_id')

        # initialize session
        Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
        session = Session()
        try:
            # delete geoserver
            try:
                geoserver = session.query(Geoserver).get(geoserver_id)
            except ObjectDeletedError:
                session.close()
                return JsonResponse({'error': "The geoserver to delete does not exist."})
            session.delete(geoserver)
            session.commit()
            session.close()
        except IntegrityError:
            session.close()
            return JsonResponse(
                {'error': "This geoserver is connected with a watershed! Must remove connection to delete."})
        return JsonResponse({'success': "Geoserver sucessfully deleted!"})
    return JsonResponse({'error': "A problem with your request exists."})

@user_passes_test(user_permission_test)
def region_delete(request):
    """
    Controller for deleting a region.
    """
    if request.is_ajax() and request.method == 'POST':
        # get/check information from AJAX request
        post_info = request.POST
        region_id = post_info.get('region_id')

        # initialize session
        Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
        session = Session()
        try:
            # delete region
            try:
                region = session.query(Region).get(region_id)
            except ObjectDeletedError:
                session.close()
                return JsonResponse({'error': "The geoserver to delete does not exist."})
            display_name = region.display_name
            region_store = ''.join(display_name.split()).lower()
            geoserver_id = region.geoserver_id
            geoserver = session.query(Geoserver).get(geoserver_id)
            geoserver_url = geoserver.url
            uname = geoserver.username
            pwd = geoserver.password

            spatial_dataset_engine = GeoServerSpatialDatasetEngine(endpoint=geoserver_url, username=uname,
                                                                   password=pwd)

            stores = spatial_dataset_engine.list_stores()

            for store in stores['result']:
                if store.endswith(region_store):
                    spatial_dataset_engine.delete_store(store,purge=True,recurse=True)

            FILE_DIR = os.path.join(GRACE_NETCDF_DIR, '')

            region_dir = os.path.join(FILE_DIR + region_store, '')

            session.delete(region)
            session.commit()
            session.close()
        except IntegrityError:
            session.close()
            return JsonResponse(
                {'error': "This geoserver is connected with a watershed! Must remove connection to delete."})
        finally:
        # Delete the temporary directory once the geojson string is created
            if region_dir is not None:
                if os.path.exists(region_dir):
                    shutil.rmtree(region_dir)
        return JsonResponse({'success': "Region sucessfully deleted!"})
    return JsonResponse({'error': "A problem with your request exists."})

def upload_shp(request):

    return_obj = {
        'success':False
    }
    GLOBAL_DIR = os.path.join(GLOBAL_NETCDF_DIR, '')

    for file in os.listdir(GLOBAL_DIR):
        if file.startswith('GRC') and file.endswith('.nc'):
            gbyos_grc_ncf = GLOBAL_NETCDF_DIR + file
        if file.startswith('CLM4') and file.endswith('.nc'):
            gbyos_fct_ncf = GLOBAL_NETCDF_DIR + file


    #Check if its an ajax post request
    if request.is_ajax() and request.method == 'POST':
        #Gettings the file list and converting the files to a geojson object. See utilities.py for the convert_shp function.
        # shapefile = request.FILES.getlist('files')
        # vals_from_shp(shapefile,GLOBAL_DIR)
        file_list = request.FILES.getlist('files')

        shp_json = convert_shp(file_list)
        gjson_obj = json.loads(shp_json)




        geometry = gjson_obj["features"][0]["geometry"]
        shape_obj = shapely.geometry.asShape(geometry)
        poly_bounds = shape_obj.bounds

        # reproj_poly_bounds = convert_shp_bounds(poly_bounds)

        #Returning the bounds and the geo_json object as a json object
        return_obj["bounds"] = poly_bounds
        return_obj["geo_json"] = gjson_obj
        return_obj["success"] = True

    return JsonResponse(return_obj)