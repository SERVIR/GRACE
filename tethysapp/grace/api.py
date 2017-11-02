from django.http import JsonResponse
from utilities import *
import json

def api_get_point_values(request):
    json_obj = {}

    if request.method == 'GET' and 'start_date' in request.GET:

        latitude = None
        longitude = None
        start_date = None
        end_date = None

        if request.GET.get('latitude'):
            latitude = request.GET['latitude']
        if request.GET.get('longitude'):
            longitude = request.GET['longitude']


        if request.GET.get('start_date'):
            start_date = request.GET['start_date']
        if request.GET.get('end_date'):
            end_date = request.GET['end_date']

        coords = str(longitude) + ',' + str(latitude)

        try:
            graph = get_global_plot_api(coords,start_date,end_date)
            graph = json.loads(graph)
            json_obj = graph

            return JsonResponse(json_obj)  # Return the json object with a list of the time and corresponding values

        except:
            json_obj = {
                "Error": "Error Processing Request"}  # Show an error if there are any problems executing the script.

            return JsonResponse(json_obj)
    else:

        latitude = None
        longitude = None

        if request.GET.get('latitude'):
            latitude = request.GET['latitude']
        if request.GET.get('longitude'):
            longitude = request.GET['longitude']


        coords = str(longitude) + ',' + str(latitude)

        try:
            graph = get_global_plot(coords)
            graph = json.loads(graph)
            json_obj = graph

            return JsonResponse(json_obj)  # Return the json object with a list of the time and corresponding values

        except:
            json_obj = {
                "Error": "Error Processing Request"}  # Show an error if there are any problems executing the script.

            return JsonResponse(json_obj)
