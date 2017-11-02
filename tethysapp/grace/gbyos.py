import sys
import os.path
import subprocess
from netCDF4 import Dataset
import netCDF4
import numpy
from datetime import datetime,timedelta
import fiona
import shapely.geometry
import math
import rtree
import csv
import tempfile, shutil
from grace import *
from .app import Grace

from django.http import JsonResponse, HttpResponse, Http404

def process_shapefile(shapefile,url,uname,pwd,region_name,GRACE_NETCDF_DIR,GLOBAL_NETCDF_DIR,display_name,geoserver_id):

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
                gbyos_grc_ncf = GLOBAL_NETCDF_DIR+file
            if file.startswith('CLM4') and file.endswith('.nc'):
                gbyos_fct_ncf = GLOBAL_NETCDF_DIR+file



        GRACE_NETCDF_DIR = os.path.join(GRACE_NETCDF_DIR,'')
        output_dir = GRACE_NETCDF_DIR+region_name

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            output_dir = os.path.join(output_dir,'')



        gbyos_wsa_csv = output_dir+region_name+".csv"
        gbyos_wsa_ncf = output_dir+region_name+".nc"

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
        gbyos_dat_str = datetime.strptime('2002-01-01T00:00:00','%Y-%m-%dT%H:%M:%S')

        YV_grc_time = []
        for JS_grc_time in range(IS_grc_time):
            gbyos_dat_dlt = timedelta(days=ZV_grc_time[JS_grc_time])
            YS_grc_time = (gbyos_dat_str + gbyos_dat_dlt).strftime('%m/%d/%Y')
            YV_grc_time.append(YS_grc_time)

        # *******************************************************************************
        # Write gbyos_wsa_csv
        # *******************************************************************************
        print('Write gbyos_wsa_csv')
        with open(gbyos_wsa_csv, 'wb') as csvfile:
            # csvwriter = csv.writer(csvfile, dialect='excel', quotechar="'",           \
            #                       quoting=csv.QUOTE_NONNUMERIC)
            csvwriter = csv.writer(csvfile, dialect='excel')
            for JS_grc_time in range(IS_grc_time):
                IV_line = [YV_grc_time[JS_grc_time], ZV_wsa[JS_grc_time]]
                csvwriter.writerow(IV_line)


        # *******************************************************************************
        # Write gbyos_wsa_ncf
        # *******************************************************************************
        print('Write gbyos_wsa_ncf')
        # -------------------------------------------------------------------------------
        # Create netCDF file
        # -------------------------------------------------------------------------------

        print('- Create netCDF file')
        h = netCDF4.Dataset(gbyos_wsa_ncf, 'w', format="NETCDF3_CLASSIC")

        time = h.createDimension("time", None)
        lat = h.createDimension("lat", IS_grc_lat)
        lon = h.createDimension("lon", IS_grc_lon)
        nv = h.createDimension("nv", 2)

        time = h.createVariable("time", "i4", ("time",))
        time_bnds = h.createVariable("time_bnds", "i4", ("time", "nv",))
        lat = h.createVariable("lat", "f4", ("lat",))
        lon = h.createVariable("lon", "f4", ("lon",))
        lwe_thickness = h.createVariable("lwe_thickness", "f4", ("time", "lat", "lon",), \
                                         fill_value=ZS_grc_fil)
        crs = h.createVariable("crs", "i4")

        # -------------------------------------------------------------------------------
        # Metadata in netCDF global attributes
        # -------------------------------------------------------------------------------
        print('- Populate global attributes')

        dt = datetime.utcnow()
        dt = dt.replace(microsecond=0)
        # Current UTC time without the microseconds
        # vsn = subprocess.Popen('../version.sh', stdout=subprocess.PIPE).communicate()
        # vsn = vsn[0]
        # vsn = vsn.rstrip()
        # Version of GBYOS

        h.Conventions = 'CF-1.6'
        h.title = ''
        h.institution = ''
        h.source = 'GBYOS: ' + ', GRACE: ' + os.path.basename(gbyos_grc_ncf) \
                   + ', Scale factors: ' + os.path.basename(gbyos_fct_ncf)
        h.history = 'date created: ' + dt.isoformat() + '+00:00'
        h.references = 'https://github.com/c-h-david/gbyos/'
        h.comment = ''
        h.featureType = 'timeSeries'

        # -------------------------------------------------------------------------------
        # Metadata in netCDF variable attributes
        # -------------------------------------------------------------------------------
        print('- Copy existing variable attributes')

        if 'time' in f.variables:
            var = f.variables['time']
            if 'standard_name' in var.ncattrs(): time.standard_name = var.standard_name
            if 'long_name' in var.ncattrs(): time.long_name = var.long_name
            if 'units' in var.ncattrs(): time.units = var.units
            if 'axis' in var.ncattrs(): time.axis = var.axis
            if 'calendar' in var.ncattrs(): time.calendar = var.calendar
            if 'bounds' in var.ncattrs(): time.bounds = var.bounds

        if 'lat' in f.variables:
            var = f.variables['lat']
            if 'standard_name' in var.ncattrs(): lat.standard_name = var.standard_name
            if 'long_name' in var.ncattrs(): lat.long_name = var.long_name
            if 'units' in var.ncattrs(): lat.units = var.units
            if 'axis' in var.ncattrs(): lat.axis = var.axis

        if 'lon' in f.variables:
            var = f.variables['lon']
            if 'standard_name' in var.ncattrs(): lon.standard_name = var.standard_name
            if 'long_name' in var.ncattrs(): lon.long_name = var.long_name
            if 'units' in var.ncattrs(): lon.units = var.units
            if 'axis' in var.ncattrs(): lon.axis = var.axis

        if 'lwe_thickness' in f.variables:
            var = f.variables['lwe_thickness']
            if 'standard_name' in var.ncattrs(): lwe_thickness.standard_name = var.standard_name
            if 'long_name' in var.ncattrs(): lwe_thickness.long_name = var.long_name
            if 'units' in var.ncattrs(): lwe_thickness.units = var.units
            if 'units' in var.ncattrs(): lwe_thickness.coordinates = var.coordinates
            if 'grid_mapping' in var.ncattrs(): lwe_thickness.grid_mapping = var.grid_mapping
            if 'cell_methods' in var.ncattrs(): lwe_thickness.cell_methods = var.cell_methods

        if 'crs' in f.variables:
            var = f.variables['crs']
            if 'grid_mapping_name' in var.ncattrs(): crs.grid_mapping_name = var.grid_mapping_name
            if 'semi_major_axis' in var.ncattrs(): crs.semi_major_axis = var.semi_major_axis
            if 'inverse_flattening' in var.ncattrs(): crs.inverse_flattening = var.inverse_flattening
        print('- Modify CRS variable attributes')
        lwe_thickness.grid_mapping = 'crs'
        crs.grid_mapping_name = 'latitude_longitude'
        crs.semi_major_axis = '6378137'
        crs.inverse_flattening = '298.257223563'
        # These are for the WGS84 spheroid

        # -------------------------------------------------------------------------------
        # Populate static data
        # -------------------------------------------------------------------------------
        print('- Populate static data')

        lon[:] = ZV_grc_lon[:]
        lat[:] = ZV_grc_lat[:]
        # Coordinates

        # -------------------------------------------------------------------------------
        # Populate dynamic data
        # -------------------------------------------------------------------------------
        print('- Populate dynamic data')

        for JS_dom_tot in range(IS_dom_tot):
            JS_grc_lon = IV_dom_lon[JS_dom_tot]
            JS_grc_lat = IV_dom_lat[JS_dom_tot]
            ZS_dom_avg = ZV_dom_avg[JS_dom_tot]
            if (ZM_grc_scl.mask[JS_grc_lat, JS_grc_lon]):
                ZS_dom_scl = 0
            else:
                ZS_dom_scl = ZM_grc_scl[JS_grc_lat, JS_grc_lon]
            for JS_grc_time in range(IS_grc_time):
                lwe_thickness[JS_grc_time, JS_grc_lat, JS_grc_lon] = \
                    f.variables['lwe_thickness'][JS_grc_time, JS_grc_lat, JS_grc_lon] \
                    - ZS_dom_avg

        time[:] = f.variables['time'][:]


        # *******************************************************************************
        # Close netCDF files
        # *******************************************************************************
        print('Close netCDF files')
        f.close()
        g.close()
        h.close()

        geotiff_output_dir = output_dir + 'geotiff'

        if not os.path.exists(geotiff_output_dir):
            os.makedirs(geotiff_output_dir)
            geotiff_output_dir = os.path.join(geotiff_output_dir, '')

        create_geotiffs(output_dir, geotiff_output_dir)
        upload_tiff(geotiff_output_dir,region_name,url,"grace",uname,pwd)
        Session = Grace.get_persistent_store_database('main_db', as_sessionmaker=True)
        session = Session()
        region = Region(geoserver_id=geoserver_id,display_name=display_name, latlon_bbox=str(bbox_val))
        session.add(region)
        session.commit()
        session.close()
        return JsonResponse({"success": "success"})
    except Exception as e:
        # Delete the temporary directory once the geojson string is created
        if output_dir is not None:
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)
        if temp_dir is not None:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
        return JsonResponse({"error":e})

    finally:
        # Delete the temporary directory once the geojson string is created
        if temp_dir is not None:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)