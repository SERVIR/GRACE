Developers API
========================

*A Developers API is provided for those who wish to incorporate the GRACE Tethys APP data into their own separate application or script.*

API Methods
-------------------

All API methods must be called using the following pattern:

``{{ host }}/apps/[parent-app]/api/[MethodName]/?param1=value1&#38;param2=value2&#38;...paramN=valueN``


GetPointValues
~~~~~~~~~~~~~~~~~~~~

Parent App: grace

Supported Methods: GET

Returns A JSON object with a timeseries for a given point

Parameters:

    latitude: Latitude in WGS 84 projection [-45,45] (Required)

    longitude: Longitude in WGS 84 projection [0,360] (Required)

    start_date: Started Date for the data [YYYY-MM-DD] (Optional)

    end_date: End Date for the data [YYYY-MM-DD] (Optional)

Example:

``http://tethys.servirglobal.net/apps/grace/api/GetPointValues/?latitude=20.7&longitude=80.2``


