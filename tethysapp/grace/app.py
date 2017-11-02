from tethys_sdk.base import TethysAppBase, url_map_maker
from tethys_sdk.app_settings import PersistentStoreDatabaseSetting, PersistentStoreConnectionSetting

class Grace(TethysAppBase):
    """
    Tethys app class for GRACE.
    """

    name = 'GRACE'
    index = 'grace:home'
    icon = 'grace/images/logo.jpg'
    package = 'grace'
    root_url = 'grace'
    color = '#e74c3c'
    description = 'View GRACE mission data'
    tags = 'Hydrology'
    enable_feedback = False
    feedback_emails = []

    def url_maps(self):
        """
        Add controllers
        """
        UrlMap = url_map_maker(self.root_url)

        url_maps = (
            UrlMap(
                name='home',
                url='grace',
                controller='grace.controllers.home'
            ),
            UrlMap(name='api',
                   url='grace/api',
                   controller='grace.controllers.api'),
            UrlMap(name='global-map',
                   url='grace/global-map',
                   controller='grace.controllers.global_map'),
            UrlMap(name='grace/plot',
                   url='grace/plot',
                   controller='grace.ajax_controllers.get_plot'),
            UrlMap(name='grace/plot-global',
                   url='grace/plot-global',
                   controller='grace.ajax_controllers.get_plot_global'),
            UrlMap(name='add-region',
                   url='grace/add-region',
                   controller='grace.controllers.add_region'),
            UrlMap(name='add-region-ajax',
                   url='grace/add-region/submit',
                   controller='grace.ajax_controllers.region_add'),
            UrlMap(name='manage-regions',
                   url='grace/manage-regions',
                   controller='grace.controllers.manage_regions'),
            UrlMap(name='manage-regions-table',
                   url='grace/manage-regions/table',
                   controller='grace.controllers.manage_regions_table'),
            UrlMap(name='delete-regions-ajax',
                   url='grace/manage-regions/delete',
                   controller='grace.ajax_controllers.region_delete'),
            UrlMap(name='add-geoserver',
                   url='grace/add-geoserver',
                   controller='grace.controllers.add_geoserver'),
            UrlMap(name='add-geoserver-ajax',
                   url='grace/add-geoserver/submit',
                   controller='grace.ajax_controllers.geoserver_add'),
            UrlMap(name='manage-geoservers',
                   url='grace/manage-geoservers',
                   controller='grace.controllers.manage_geoservers'),
            UrlMap(name='manage-geoservers-table',
                   url='grace/manage-geoservers/table',
                   controller='grace.controllers.manage_geoservers_table'),
            UrlMap(name='update-geoservers-ajax',
                   url='grace/manage-geoservers/submit',
                   controller='grace.ajax_controllers.geoserver_update'),
            UrlMap(name='delete-geoserver-ajax',
                   url='grace/manage-geoservers/delete',
                   controller='grace.ajax_controllers.geoserver_delete'),
            UrlMap(name='map',
                   url='grace/map',
                   controller='grace.controllers.map'),
            UrlMap(name='plot-region',
                   url='grace/plot-region',
                   controller='grace.ajax_controllers.plot_region'),
            UrlMap(name='api_get_point_values',
                   url='grace/api/GetPointValues',
                   controller='grace.api.api_get_point_values'),
            UrlMap(name='upload-shp',
                   url='grace/upload-shp',
                   controller='grace.ajax_controllers.upload_shp'),
        )

        return url_maps


    def persistent_store_settings(self):

        ps_settings = (
            PersistentStoreDatabaseSetting(
                name='main_db',
                description='For storing Region and Geoserver metadata',
                required=True,
                initializer='grace.model.init_main_db',
                spatial=False,
            ),
        )

        return ps_settings
