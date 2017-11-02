Web Application
==================

**This app is created to run in the Teyths Platform programming environment.** 


**You can find a working demo here at ** http://tethys.servirglobal.net/apps/grace/

Prerequisites
--------------

- Tethys Platform (CKAN, PostgresQL, GeoServer) with GeoServer
- gdal (Python package for working with geospatial data)
- numpy (Python package for scientific computing)
- netCDF4 (Python package for reading and writing netCDF files)
- fiona (Python package for reading and writing spatial data files)
- shapely (Python package for set-theoretic analysis and manipulation of planar features)
- pyshp (Python package for reading and writing ESRI shapefiles)
- rtree (Python wrapper of libpspatialindex)
- GeoServer needs CORS enabled
- GRACE Tellus and CLM4 NetCDF files. See: ftp://podaac.jpl.nasa.gov/allData/tellus/L3/mascon/RL05/JPL/CRI/netcdf/
 
Install Tethys Platform
~~~~~~~~~~~~~~~~~~~~~~~~~~~
See: http://docs.tethysplatform.org/en/latest/installation.html


Install the python packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note: Before installing the python packages into your python site-packages, activate
your Tethys conda environment using the alias `t`:

::

    $ t

::

    (tethys)$ conda install -c conda-forge gdal numpy netCDF4 fiona shapely pyshp rtree


Download GRACE Global Datasets
---------------------------------
Prior to installing the app, download the GRACE NetCDF files from the JPL FTP site at ftp://podaac.jpl.nasa.gov/allData/tellus/L3/mascon/RL05/JPL/CRI/netcdf/

.. _ftp://podaac.jpl.nasa.gov/allData/tellus/L3/mascon/RL05/JPL/CRI/netcdf/: ftp://podaac.jpl.nasa.gov/allData/tellus/L3/mascon/RL05/JPL/CRI/netcdf/

Download the following files only:
    *CLM4.SCALE_FACTOR.JPL.MSCNv01CRIv01.nc*
    *GRCTellus.JPL.200204_201706.GLO.RL05M_1.MSCNv02CRIv02.nc*

Place them in a direcotry called grace with subdirectory named global as follows.

::

    $ cd /
    $ mkdir grace
    $ cd grace
    $ tethys@tethys:/grace $ cd global
    $ tethys@tethys:/grace/global wget ftp://podaac.jpl.nasa.gov/allData/tellus/L3/mascon/RL05/JPL/CRI/netcdf/GRCTellus.JPL.200204_201706.GLO.RL05M_1.MSCNv02CRIv02.nc
    $ tethys@tethys:/grace/global wget ftp://podaac.jpl.nasa.gov/allData/tellus/L3/mascon/RL05/JPL/CRI/netcdf/CLM4.SCALE_FACTOR.JPL.MSCNv01CRIv01.nc

.. NOTE::
    The GRCTellus file is updated monthly. As of now there isn't a way to automatically update the GRACE Tellus file within the app. This feature will be added in a future release. But for now you will need to download this file using curl or wget onto your local/production machine.

Change the permissions on the grace folder to be readable, writable and executable


::

    $ sudo chmod -R 777 /home/mymachine/grace/

.. WARNING::
    Not changing the permissions will cause the app's subsetting functionality to fail

Web App Installation
----------------------

Installation for App Development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Download the source code from github

::

    $ t
    (tethys)$ git clone https://github.com/BYU-Hydroinformatics/grace
    (tethys)$ cd grace
    (tethys)$ python setup.py develop

Change the NetCDF file paths in controllers.py and ajax_controllers.py to the path of the GRACE directory in your local machine

Open the :file:`controllers.py` for editing using ``vim`` or any text editor of your choice:

::

    (tethys)$ cd grace/tethysapp/grace
    (tethys)$ sudo vi controllers.py

Press :kbd:`i` to start editing and change the ``NETCDF_DIR`` global variable. You can find it right after the import statements.

::

    NETCDF_DIR = '/home/mymachine/grace/'

Press :kbd:`ESC` to exit ``INSERT`` mode and then press ``:x`` and :kbd:`ENTER` to save changes and exit.

Open the :file:`ajax_controllers.py` for editing using ``vim`` or any text editor of your choice:

::

    (tethys)$ sudo vi controllers.py

Press :kbd:`i` to start editing and change the ``GRACE_NETCDF_DIR`` and ``GLOBAL_NETCDF_DIR`` global variables. You can find them right after the import statements.

::

    GRACE_NETCDF_DIR = '/home/mymachine/grace/'
    GLOBAL_NETCDF_DIR = '/home/mymachine/grace/global/'

Press :kbd:`ESC` to exit ``INSERT`` mode and then press ``:x`` and :kbd:`ENTER` to save changes and exit.

Start the Tethys Server

::

    (tethys)$ tms

You should now have the GRACE Viewer running on a development server on your machine. Tethys Platform provides a web interface called the Tethys Portal. You can access the app through the Tethys portal by opening http://localhost:8000/ (or if you provided custom host and port options to the install script then it will be <HOST>:<PORT>) in a new tab in your web browser.

Installation for Production
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Installing apps in a Tethys Platform configured for production can be challenging. Most of the difficulties arise, because Tethys is served by Nginx in production and all the files need to be owned by the Nginx user. The following instructions will allow you to deploy the GRACE Viewer on your own Tethys production server. You can find the Tethys Production installation instructions `here. <http://docs.tethysplatform.org/en/stable/installation/production.html>`_

Change the Ownership of the Files to the Current User

*During the production installation any Tethys related files were change to be owned by the Nginx user. To make any changes on the server it is easiest to change the ownership back to the current user. This is easily done with an alias that was created in the tethys environment during the production installation process*


::

    $ t
    (tethys)$ tethys_user_own

Download App Source Code from GitHub

::

    $ cd $TETHYS_HOME/apps/
    $ sudo git clone https://github.com/BYU-Hydroinformatics/grace

.. tip::

    Substitute $TETHYS_HOME with the path to the tethys main directory.


Change the NetCDF file paths in controllers.py and ajax_controllers.py to the path of the GRACE directory in your local machine

Open the :file:`controllers.py` for editing using ``vim`` or any text editor of your choice:

::

    (tethys)$ cd $TETHYS_HOME/apps/grace/tethysapp/grace
    (tethys)$ sudo vi controllers.py

Press :kbd:`i` to start editing and change the ``NETCDF_DIR`` global variable. You can find it right after the import statements.

::

    NETCDF_DIR = 'home/mymachine/grace/'

Press :kbd:`ESC` to exit ``INSERT`` mode and then press ``:x`` and :kbd:`ENTER` to save changes and exit.

Open the :file:`ajax_controllers.py` for editing using ``vim`` or any text editor of your choice:

::

    (tethys)$ sudo vi controllers.py

Press :kbd:`i` to start editing and change the ``GRACE_NETCDF_DIR`` and ``GLOBAL_NETCDF_DIR`` global variables. You can find them right after the import statements.

::

    GRACE_NETCDF_DIR = '/home/mymachine/grace/'
    GLOBAL_NETCDF_DIR = '/home/mymachine/grace/global/'

Press :kbd:`ESC` to exit ``INSERT`` mode and then press ``:x`` and :kbd:`ENTER` to save changes and exit.

Return to the main directory of the app. Then, execute the setup script (:file:`setup.py`) with the ``install`` command to make Python aware of the app and install any of its dependencies

::

    (tethys)$ cd $TETHYS_HOME/apps/grace/
    (tethys)$ python setup.py install

Collect Static Files and Workspaces

The static files and files in app workspaces are hosted by Nginx, which necessitates collecting all of the static files to a single directory and all workspaces to another single directory. These directory is configured through the ``STATIC_ROOT`` and ``TETHYS_WORKSPACES_ROOT`` setting in the :file:`settings.py` file. Collect the static files and workspaces with this command

::

    (tethys)$ tethys manage collectall


Change the Ownership of Files to the Nginx User

The Nginx user must own any files that Nginx is serving. This includes the source files, static files, and any workspaces that your app may have. The following alias will accomplish the change in ownership that is required

::

    (tethys)$ tethys_server_own
     

Restart uWSGI and Nginx services to effect the changes

::

    $ sudo systemctl restart tethys.uwsgi.service
    $ sudo systemctl restart nginx

.. note::

    For updating the app on production server, simply pull the app from GitHub. Once you have made a pull request (at times you may have to stash your local changes), follow the above steps to reinstall/update the app.


Enable CORS on geoserver
--------------------------
The following workflow is necessary if the geoserver is not on the same server as the Tethys platform.


Create a new bash session in the tethys_geoserver docker container:

::

    $ t
    (tethys)$ docker exec -it tethys_geoserver /bin/bash
    (tethys)$ cd node1/webapps/geoserver/WEB-INF


.. NOTE::
    You can make this change to any other node in the geoserver docker.


Insert the following in the filters list:

::

    <filter>
    <filter-name>CorsFilter</filter-name>
    <filter-class>org.apache.catalina.filters.CorsFilter</filter-class>
    <init-param>
      <param-name>cors.allowed.origins</param-name>
      <param-value>http://127.0.0.1:8000, http://127.0.0.1:8181</param-value>
    </init-param>
    </filter>


Insert this filter-mapping to the filter-mapping list:

::

    <filter-mapping>
    <filter-name>CorsFilter</filter-name>
    <url-pattern>/*</url-pattern>
    </filter-mapping>


Save the web.xml file

::

    $ exit
    $ docker restart tethys_geoserver
