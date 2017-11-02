from django.shortcuts import render
from django.contrib.auth.decorators import login_required,user_passes_test
from tethys_sdk.gizmos import *
import csv, os
from datetime import datetime,timedelta
from tethys_sdk.services import get_spatial_dataset_engine
import urlparse
from grace import *
from utilities import *
import json,time
from .app import Grace
from model import *

NETCDF_DIR = '/grace/'


def home(request):
    """
    Controller for the app home page.
    """

    #create_global_tiff("/grace/global/GRCTellus.JPL.200204_201608.GLO.RL05M_1.MSCNv02CRIv02.nc","/home/tethys/geotiff_global/","lwe_thickness")
    #upload_global_tiff("/home/tethys/geotiff_global/","http://tethys.servirglobal.net:8181/geoserver/rest/","graceglobal")
    Session = Grace.get_persistent_store_database('main_db',as_sessionmaker=True)
    session = Session()
    # Query DB for regions
    regions = session.query(Region).all()
    region_list = []

    for region in regions:
        region_list.append(("%s" % (region.display_name), region.id))

    session.close()
    if region_list:
        region_select = SelectInput(display_text='Select a Region',
                                    name='region-select',
                                    options=region_list, )
    else:
        region_select = None

    context = {
        "region_select": region_select, "regions_length": len(region_list), 'host': 'http://%s' % request.get_host()
    }

    return render(request, 'grace/home.html', context)


def api(request):

    context = {'host': 'http://%s' % request.get_host()}

    return render(request, 'grace/api.html', context)




def map(request):

    context = {}

    info = request.GET

    region_id = info.get('region-select')
    Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
    session = Session()

    region = session.query(Region).get(region_id)
    display_name = region.display_name

    bbox = [float(x) for x in region.latlon_bbox.strip("(").strip(")").split(',')]
    json.dumps(bbox)

    geoserver = session.query(Geoserver).get(region.geoserver_id)
    geoserver_url = geoserver.url
    region_store = ''.join(display_name.split()).lower()


    FILE_DIR = os.path.join(NETCDF_DIR,'')

    region_dir = os.path.join(FILE_DIR+region_store,'')

    geotiff_dir = os.path.join(region_dir+"geotiff")

    sorted_files = sorted(os.listdir(geotiff_dir), key=lambda x: datetime.strptime(x, '%Y_%m_%d.tif'))
    layers_length = len(sorted_files)
    grace_layer_options = []

    for file in sorted_files:
        year = int(file[:-4].split('_')[0])
        month = int(file[:-4].split('_')[1])
        day = int(file[:-4].split('_')[2])
        date_str = datetime(year,month,day)
        date_str = date_str.strftime("%Y %B %d")
        grace_layer_options.append([date_str,file[:-4]+"_"+region_store])

    select_layer = SelectInput(display_text='Select a day',
                               name='select_layer',
                               multiple=False,
                               options=grace_layer_options, )

    csv_file = region_dir+region_store+".csv"
    with open(csv_file, 'rb') as f:
        reader = csv.reader(f)
        csvlist = list(reader)

    volume_time_series = []
    volume = []
    x_tracker = []
    formatter_string = "%m/%d/%Y"
    for item in csvlist:
        mydate = datetime.strptime(item[0], formatter_string)
        mydate = time.mktime(mydate.timetuple()) * 1000
        volume_time_series.append([mydate, float(item[1])])
        volume.append(float(item[1]))
        x_tracker.append(mydate)

    range = [round(min(volume), 2), round(max(volume), 2)]
    range = json.dumps(range)

    # Configure the time series Plot View
    grace_plot = TimeSeries(
        engine='highcharts',
        title=display_name+ ' GRACE Data',
        y_axis_title='Total Terrestrial Water Storage Anomaly',
        y_axis_units='cm',
        series=[
            {
                'name': 'Height of Liquid Water',
                'color': '#0066ff',
                'data': volume_time_series,
            },
            {
                'name': 'Tracker',
                'color': '#ff0000',
                'data': [[min(x_tracker), round(min(volume), 2)], [min(x_tracker), round(max(volume), 2)]]
            },
        ],
        width='100%',
        height='300px'
    )

    wms_url = geoserver_url[:-5]+"wms"
    color_bar = get_color_bar()
    color_bar = json.dumps(color_bar)

    if bbox[0] < 0 and bbox[2] < 0:
        map_center = [( (360+(int(bbox[0])))+(360+(int(bbox[2])))) / 2,(int(bbox[1])+int(bbox[3])) / 2]
    else:
        map_center = [(int(bbox[0]) + int(bbox[2])) / 2, (int(bbox[1]) + int(bbox[3])) / 2]
    json.dumps(map_center)
    json.dumps(x_tracker)

    context = {"region_id":region_id,"display_name":display_name,"wms_url":wms_url,"select_layer":select_layer,"layers_length":layers_length,"grace_plot":grace_plot
               ,'x_tracker':x_tracker,"color_bar":color_bar,"range":range,"bbox":bbox,"map_center":map_center}

    return render(request, 'grace/map.html', context)

# @login_required
# def nepal_graph(request):
#
#     #Creating the Chart
#
#     user_workspace = Grace.get_app_workspace()
#
#     csv_file = os.path.join(user_workspace.path, 'output/Nepal/hydrograph.csv')
#
#
#
#     with open(csv_file, 'rb') as f:
#         reader = csv.reader(f)
#         csvlist = list(reader)
#
#     volume_time_series = []
#     volume = []
#     x_tracker = []
#     formatter_string = "%m/%d/%Y"
#     for item in csvlist:
#         mydate = datetime.strptime(item[0], formatter_string)
#         mydate = time.mktime(mydate.timetuple())*1000
#         volume_time_series.append([mydate, float(item[1])])
#         volume.append(float(item[1]))
#         x_tracker.append(mydate)
#
#     range = [round(min(volume),2),round(max(volume),2)]
#     range = json.dumps(range)
#
#     # Configure the time series Plot View
#     grace_plot = TimeSeries(
#         engine='highcharts',
#         title= 'Nepal GRACE Data',
#         y_axis_title='Volume',
#         y_axis_units='cm',
#         series=[
#             {
#                 'name': 'Change in Volume',
#                 'color': '#0066ff',
#                 'data': volume_time_series,
#             },
#             {
#                 'name':'Tracker',
#                 'color': '#ff0000',
#                 'data':[[min(x_tracker),-50],[min(x_tracker),50]]
#             },
#         ],
#         width='100%',
#         height='300px'
#     )
#
#     #Connecting to the Geoserver
#     geoserver_engine = get_spatial_dataset_engine(name='default')
#     stores = geoserver_engine.list_stores(workspace='grace')
#
#     grace_layer_options = []
#     sorted_stores = sorted(stores['result'],key=lambda x:datetime.strptime(x,'%Y_%m_%d_nepal'))
#     for store in sorted_stores:
#
#         year = int(store.split('_')[0])
#         month = int(store.split('_')[1])
#         day = int(store.split('_')[2])
#         date_str = datetime(year,month,day)
#         date_str = date_str.strftime("%Y %B %d")
#         grace_layer_options.append([date_str,store])
#
#     slider_max = len(grace_layer_options)
#
#     select_layer = SelectInput(display_text='Select a day',
#                                name='select_layer',
#                                multiple=False,
#                                options=grace_layer_options,)
#     legend_file = os.path.join(user_workspace.path, 'output/Nepal/legend.csv')
#
#
#     with open(legend_file, 'rb') as f:
#         reader = csv.reader(f)
#         legend_list = list(reader)
#
#     legend_json = json.dumps(legend_list)
#     x_tracker = json.dumps(x_tracker)
#
#     color_bar = get_color_bar()
#     color_bar = json.dumps(color_bar)
#
#
#     context = {'grace_plot': grace_plot,'select_layer':select_layer,'layers_json':legend_json,'range':range,'slider_max':slider_max,'x_tracker':x_tracker,"color_bar":color_bar}
#     return render(request, 'grace/nepal_graph.html', context)


def global_map(request):

    color_bar = get_color_bar()
    color_bar = json.dumps(color_bar)
    # Connecting to the Geoserver
    # geoserver_engine = get_spatial_dataset_engine(name='grace')
    # stores = geoserver_engine.list_stores(workspace='grace')
    #
    # grace_layer_options = []
    # sorted_stores = sorted(stores['result'], key=lambda x: datetime.strptime(x, '%Y_%m_%d'))
    # for store in sorted_stores:
    #     year = int(store.split('_')[0])
    #     month = int(store.split('_')[1])
    #     day = int(store.split('_')[2])
    #     date_str = datetime(year, month, day)
    #     date_str = date_str.strftime("%Y %B %d")
    #     grace_layer_options.append([date_str, store])

    FILE_DIR = os.path.join(NETCDF_DIR, '')

    region_dir = os.path.join(FILE_DIR + 'niger', '')

    geotiff_dir = os.path.join(region_dir + "geotiff")

    sorted_files = sorted(os.listdir(geotiff_dir), key=lambda x: datetime.strptime(x, '%Y_%m_%d.tif'))
    layers_length = len(sorted_files)
    grace_layer_options = []

    for file in sorted_files:
        year = int(file[:-4].split('_')[0])
        month = int(file[:-4].split('_')[1])
        day = int(file[:-4].split('_')[2])
        date_str = datetime(year, month, day)
        date_str = date_str.strftime("%Y %B %d")
        grace_layer_options.append([date_str, file[:-4]])



    slider_max = len(grace_layer_options)

    select_layer = SelectInput(display_text='Select a day',
                               name='select_layer',
                               multiple=False,
                               options=grace_layer_options, )



    context = {'select_layer':select_layer,'slider_max':slider_max,"color_bar":color_bar}

    return render(request, 'grace/global_map.html', context)

@user_passes_test(user_permission_test)
def add_region(request):

    region_name_input = TextInput(display_text='Region Display Name',
                                     name='region-name-input',
                                     placeholder='e.g.: Utah',
                                     icon_append='glyphicon glyphicon-home',
                                     ) #Input for the Region Display Name

    Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
    session = Session()
    # Query DB for geoservers
    geoservers = session.query(Geoserver).all()
    geoserver_list = []
    for geoserver in geoservers:
        geoserver_list.append(( "%s (%s)" % (geoserver.name, geoserver.url),
                               geoserver.id))

    session.close()
    if geoserver_list:
        geoserver_select = SelectInput(display_text='Select a Geoserver',
                                       name='geoserver-select',
                                       options=geoserver_list,)
    else:
        geoserver_select = None

    add_button = Button(display_text='Add Region',
                        icon='glyphicon glyphicon-plus',
                        style='success',
                        name='submit-add-region',
                        attributes={'id': 'submit-add-region'}, )  # Add region button

    context = {"region_name_input":region_name_input, "geoserver_select": geoserver_select,"add_button":add_button}
    return render(request, 'grace/add_region.html', context)

@user_passes_test(user_permission_test)
def add_geoserver(request):
    """
        Controller for the app add_geoserver page.
    """

    geoserver_name_input = TextInput(display_text='Geoserver Name',
                                     name='geoserver-name-input',
                                     placeholder='e.g.: BYU Geoserver',
                                     icon_append='glyphicon glyphicon-tag', )

    geoserver_url_input = TextInput(display_text='Geoserver REST Url',
                                    name='geoserver-url-input',
                                    placeholder='e.g.: http://tethys.byu.edu:8181/geoserver/rest',
                                    icon_append='glyphicon glyphicon-cloud-download')

    geoserver_username_input = TextInput(display_text='Geoserver Username',
                                         name='geoserver-username-input',
                                         placeholder='e.g.: admin',
                                         icon_append='glyphicon glyphicon-user', )

    add_button = Button(display_text='Add Geoserver',
                        icon='glyphicon glyphicon-plus',
                        style='success',
                        name='submit-add-geoserver',
                        attributes={'id': 'submit-add-geoserver'}, )

    context = {
        'geoserver_name_input': geoserver_name_input,
        'geoserver_url_input': geoserver_url_input,
        'geoserver_username_input': geoserver_username_input,
        'add_button': add_button,
    }

    return render(request, 'grace/add_geoserver.html', context)

@user_passes_test(user_permission_test)
def manage_regions(request):
    """
    Controller for the app manage_geoservers page.
    """
    #initialize session
    Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
    session = Session()
    num_regions = session.query(Region).count()

    session.close()

    context = {
                'initial_page': 0,
                'num_regions': num_regions,
              }

    return render(request, 'grace/manage_regions.html', context)

@user_passes_test(user_permission_test)
def manage_regions_table(request):
    """
    Controller for the app manage_geoservers page.
    """
    #initialize session
    Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
    session = Session()
    RESULTS_PER_PAGE = 5
    page = int(request.GET.get('page'))

    # Query DB for data store types
    regions = session.query(Region)\
                        .order_by(Region.display_name) \
                        .all()[(page * RESULTS_PER_PAGE):((page + 1)*RESULTS_PER_PAGE)]

    prev_button = Button(display_text='Previous',
                         name='prev_button',
                         attributes={'class':'nav_button'},)

    next_button = Button(display_text='Next',
                         name='next_button',
                         attributes={'class':'nav_button'},)

    context = {
                'prev_button' : prev_button,
                'next_button': next_button,
                'regions': regions,
              }

    session.close()

    return render(request, 'grace/manage_regions_table.html', context)
@user_passes_test(user_permission_test)
def manage_geoservers(request):
    """
    Controller for the app manage_geoservers page.
    """
    #initialize session
    Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
    session = Session()
    num_geoservers = session.query(Geoserver).count()
    session.close()

    context = {
                'initial_page': 0,
                'num_geoservers': num_geoservers,
              }

    return render(request, 'grace/manage_geoservers.html', context)

@user_passes_test(user_permission_test)
def manage_geoservers_table(request):
    """
    Controller for the app manage_geoservers page.
    """
    #initialize session
    Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
    session = Session()
    RESULTS_PER_PAGE = 5
    page = int(request.GET.get('page'))

    # Query DB for data store types
    geoservers = session.query(Geoserver)\
                        .order_by(Geoserver.name, Geoserver.url) \
                        .all()[(page * RESULTS_PER_PAGE):((page + 1)*RESULTS_PER_PAGE)]

    prev_button = Button(display_text='Previous',
                         name='prev_button',
                         attributes={'class':'nav_button'},)

    next_button = Button(display_text='Next',
                         name='next_button',
                         attributes={'class':'nav_button'},)

    context = {
                'prev_button' : prev_button,
                'next_button': next_button,
                'geoservers': geoservers,
              }

    session.close()

    return render(request, 'grace/manage_geoservers_table.html', context)



