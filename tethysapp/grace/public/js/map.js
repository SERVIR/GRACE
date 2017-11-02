/*****************************************************************************
 * FILE:    GRACE VIEWER MAP.JS
 * DATE:    22 May 2017
 * AUTHOR: Sarva Pulla
 * COPYRIGHT: (c) Brigham Young University 2017
 * LICENSE: BSD 2-Clause
 *****************************************************************************/

/*****************************************************************************
 *                      LIBRARY WRAPPER
 *****************************************************************************/

var GRACE_MAP = (function() {
    // Wrap the library in a package function
    "use strict"; // And enable strict mode for this library

    /************************************************************************
     *                      MODULE LEVEL / GLOBAL VARIABLES
     *************************************************************************/
    var animationDelay,
        bbox,
        chart,
        cb_min,
        cb_max,
        color_bar,
        current_layer,
        element,
        $get_plot,
        layers_length,
        layers,
        layers_dict,
        map,
        map_center,
        $maxSlider,
        popup,
        plotter,
        public_interface, // Object returned by the module
        $region_element,
        range,
        range_min,
        range_max,
        sliderInterval,
        shp_source,
        shp_layer,
        tracker,
        wms_url,
        wms_source,
        wms_layer;

    /************************************************************************
     *                    PRIVATE FUNCTION DECLARATIONS
     *************************************************************************/

    var add_wms,animate,clear_coords,cbar_str,gen_color_bar,get_plot,init_events,init_map,init_vars,init_slider,update_color_bar,update_wms;
    /************************************************************************
     *                    PRIVATE FUNCTION IMPLEMENTATIONS
     *************************************************************************/

    init_vars = function(){
        $region_element= $('#region');
        bbox =  $region_element.attr('data-bbox');
        bbox = JSON.parse(bbox);
        layers_length = $region_element.attr('data-layers-length');
        color_bar = $region_element.attr('data-color-bar');
        color_bar = JSON.parse(color_bar);
        range = $region_element.attr('data-range');
        range = JSON.parse(range);
        range_max = range[1];
        range_min = range[0];
        cb_min = $region_element.attr('data-range-min');
        cb_max = $region_element.attr('data-range-max');
        map_center = $region_element.attr('data-map-center');
        map_center = JSON.parse(map_center);
        plotter = $('#plotter');
        $get_plot = $('#get-plot');
        tracker = $region_element.attr('data-tracker');
        tracker = JSON.parse(tracker);
        chart = $(".highcharts-plot").highcharts();
        wms_url = $region_element.attr('data-wms-url');
        animationDelay = 1000;
        sliderInterval = {};
    };

    gen_color_bar = function(){
        var cv  = document.getElementById('cv'),
            ctx = cv.getContext('2d');
        color_bar.forEach(function(color,i){
            ctx.beginPath();
            ctx.fillStyle = color[0];
            ctx.fillRect(i*35,0,35,20);
            ctx.fillText(color[1],i*35,30);
        });

    };

    update_color_bar = function(){
        var cv  = document.getElementById('cv'),
            ctx = cv.getContext('2d');
        ctx.clearRect(0,0,cv.width,cv.height);
        color_bar.forEach(function(color,i){
            ctx.beginPath();
            ctx.fillStyle = color[0];
            ctx.fillRect(i*35,0,35,20);
            ctx.fillText(color[1],i*35,30);
        });
    };

    clear_coords = function(){
        $("#poly-lat-lon").val('');
        $("#point-lat-lon").val('');
        $("#shp-lat-lon").val('');
    };

    init_map = function(){
        var projection = ol.proj.get('EPSG:3857');
        var baseLayer = new ol.layer.Tile({
            source: new ol.source.BingMaps({
                key: '5TC0yID7CYaqv3nVQLKe~xWVt4aXWMJq2Ed72cO4xsA~ApdeyQwHyH_btMjQS1NJ7OHKY8BK-W-EMQMrIavoQUMYXeZIQOUURnKGBOC7UCt4',
                imagerySet: 'AerialWithLabels' // Options 'Aerial', 'AerialWithLabels', 'Road'
            })
        });
        var fullScreenControl = new ol.control.FullScreen();

        var view = new ol.View({
            center: ol.proj.transform(map_center,'EPSG:4326','EPSG:3857'),
            projection: projection,
            zoom: 5
        });
        //Creating an empty source and empty layer for displaying the shpefile object
        shp_source = new ol.source.Vector();
        shp_layer = new ol.layer.Vector({
            source: shp_source
        });

        var vector_source = new ol.source.Vector({
            wrapX: false
        });

        var vector_layer = new ol.layer.Vector({
            name: 'my_vectorlayer',
            source: vector_source,
            style: new ol.style.Style({
                fill: new ol.style.Fill({
                    color: 'rgba(255, 255, 255, 0.2)'
                }),
                stroke: new ol.style.Stroke({
                    color: '#ffcc33',
                    width: 2
                }),
                image: new ol.style.Circle({
                    radius: 7,
                    fill: new ol.style.Fill({
                        color: '#ffcc33'
                    })
                })
            })
        });
        wms_source = new ol.source.TileWMS();

        wms_layer = new ol.layer.Tile({
            source: wms_source
        });

        layers = [baseLayer,vector_layer,shp_layer,wms_layer];

        layers_dict = {};


        map = new ol.Map({
            target: document.getElementById("map"),
            layers: layers,
            view: view
        });
        map.addControl(new ol.control.ZoomSlider());
        map.addControl(fullScreenControl);
        map.crossOrigin = 'anonymous';
        element = document.getElementById('popup');

        popup = new ol.Overlay({
            element: element,
            positioning: 'bottom-center',
            stopEvent: true
        });

        map.addOverlay(popup);

        //Code for adding interaction for drawing on the map
        var lastFeature, draw, featureType;

        //Clear the last feature before adding a new feature to the map
        var removeLastFeature = function () {
            if (lastFeature) vector_source.removeFeature(lastFeature);
        };

        //Add interaction to the map based on the selected interaction type
        var addInteraction = function (geomtype) {
            var typeSelect = document.getElementById('types');
            var value = typeSelect.value;
            $('#data').val('');
            if (value !== 'None') {
                if (draw)
                    map.removeInteraction(draw);

                draw = new ol.interaction.Draw({
                    source: vector_source,
                    type: geomtype
                });


                map.addInteraction(draw);
            }
            if (featureType === 'Point' || featureType === 'Polygon') {

                draw.on('drawend', function (e) {
                    lastFeature = e.feature;

                });

                draw.on('drawstart', function (e) {
                    vector_source.clear();
                });

            }

        };

        vector_layer.getSource().on('addfeature', function(event){
            //Extracting the point/polygon values from the drawn feature
            var feature_json = saveData();
            var parsed_feature = JSON.parse(feature_json);
            var feature_type = parsed_feature["features"][0]["geometry"]["type"];
            if (feature_type == 'Point'){
                var coords = parsed_feature["features"][0]["geometry"]["coordinates"];
                var proj_coords = ol.proj.transform(coords, 'EPSG:3857','EPSG:4326');
                $("#point-lat-lon").val(proj_coords);

            } else if (feature_type == 'Polygon'){
                var coords = parsed_feature["features"][0]["geometry"]["coordinates"][0];
                proj_coords = [];
                coords.forEach(function (coord) {
                    var transformed = ol.proj.transform(coord,'EPSG:3857','EPSG:4326');
                    proj_coords.push('['+transformed+']');
                });
                var json_object = '{"type":"Polygon","coordinates":[['+proj_coords+']]}';
                $("#poly-lat-lon").val(json_object);
            }
        });

        function saveData() {
            // get the format the user has chosen
            var data_type = 'GeoJSON',
                // define a format the data shall be converted to
                format = new ol.format[data_type](),
                // this will be the data in the chosen format
                data;
            try {
                // convert the data of the vector_layer into the chosen format
                data = format.writeFeatures(vector_layer.getSource().getFeatures());
            } catch (e) {
                // at time of creation there is an error in the GPX format (18.7.2014)
                $('#data').val(e.name + ": " + e.message);
                return;
            }
            // $('#data').val(JSON.stringify(data, null, 4));
            return data;

        }

        $('#types').change(function (e) {
            featureType = $(this).find('option:selected').val();
            if(featureType == 'None'){
                $('#data').val('');
                clear_coords();
                map.removeInteraction(draw);
                vector_layer.getSource().clear();
                shp_layer.getSource().clear();
            }else if(featureType == 'Upload')
            {
                clear_coords();
                vector_layer.getSource().clear();
                shp_layer.getSource().clear();
                map.removeInteraction(draw);
                $modalUpload.modal('show');
            }else if(featureType == 'Point')
            {
                clear_coords();
                shp_layer.getSource().clear();
                addInteraction(featureType);
            }else if(featureType == 'Polygon'){
                clear_coords();
                shp_layer.getSource().clear();
                addInteraction(featureType);
            }
        }).change();

    };


    init_events = function() {
        (function () {
            var target, observer, config;
            // select the target node
            target = $('#app-content-wrapper')[0];

            observer = new MutationObserver(function () {
                window.setTimeout(function () {
                    map.updateSize();
                }, 350);
            });
            $(window).on('resize', function () {
                map.updateSize();
            });

            config = {attributes: true};

            observer.observe(target, config);
        }());

        map.on("singleclick",function(evt){

            $(element).popover('destroy');


            if (map.getTargetElement().style.cursor == "pointer" && $("#types").find('option:selected').val()=="None") {
                var clickCoord = evt.coordinate;
                popup.setPosition(clickCoord);
                var view = map.getView();
                var viewResolution = view.getResolution();

                var wms_url = current_layer.getSource().getGetFeatureInfoUrl(evt.coordinate, viewResolution, view.getProjection(), {'INFO_FORMAT': 'application/json'}); //Get the wms url for the clicked point
                if (wms_url) {
                    //Retrieving the details for clicked point via the url
                    $.ajax({
                        type: "GET",
                        url: wms_url,
                        dataType: 'json',
                        success: function (result) {
                            var value = parseFloat(result["features"][0]["properties"]["GRAY_INDEX"]);
                            value = value.toFixed(2);
                            $(element).popover({
                                'placement': 'top',
                                'html': true,
                                //Dynamically Generating the popup content
                                'content':'Value: '+value
                            });

                            $(element).popover('show');
                            $(element).next().css('cursor', 'text');


                        },
                        error: function (XMLHttpRequest, textStatus, errorThrown) {
                            console.log(Error);
                        }
                    });
                }
            }
        });

        map.on('pointermove', function(evt) {
            if (evt.dragging) {
                return;
            }
            var pixel = map.getEventPixel(evt.originalEvent);
            var hit = map.forEachLayerAtPixel(pixel, function(layer) {
                if (layer != layers[0]&& layer != layers[1] && layer != layers[2] && layer != layers[4]){
                    current_layer = layer;
                    return true;}
            });
            map.getTargetElement().style.cursor = hit ? 'pointer' : '';
        });
    };

    init_slider = function(){

        $( "#slider" ).slider({
            value:1,
            min: 0,
            max: layers_length - 1,
            step: 1, //Assigning the slider step based on the depths that were retrieved in the controller
            animate:"fast",
            slide: function( event, ui ) {
                var date_text = $("#select_layer option")[ui.value].text;
                $( "#grace-date" ).val(date_text); //Get the value from the slider
                var date_value = $("#select_layer option")[ui.value].value;

                update_wms(date_value);
            }
        });

        $( "#opacity-slider" ).slider({
            value:0.7,
            min: 0.2,
            max: 1,
            step: 0.1, //Assigning the slider step based on the depths that were retrieved in the controller
            animate:"fast",
            slide: function( event, ui ) {
                var opacity = ui.value;
                $("#opacity").val(opacity);
                var date_idx = $("#slider").slider("option","value");
                var date_value = $("#select_layer option")[date_idx].value;
                color_bar.forEach(function(color,i){
                    color[2] = opacity;
                });
                update_wms(date_value);
            }
        });

        $( "#max-slider" ).slider({
            value:50,
            min: 1,
            max: 50,
            step: 1, //Assigning the slider step based on the depths that were retrieved in the controller
            animate:"fast",
            slide: function( event, ui ) {
                var color_range = ui.value;
                $("#cbar-slider").val(color_range);
                var date_idx = $("#slider").slider("option","value");
                var date_value = $("#select_layer option")[date_idx].value;
                cb_max = ui.value;

                var iter_size = cb_max / 10;
                var cbar_val = -cb_max;
                var new_cbar = [];
                for (var i=0;i<=20;i+=1){
                    new_cbar.push(parseFloat(cbar_val).toFixed(1));
                    cbar_val += iter_size;
                }
                color_bar.forEach(function(color,i){
                    color[1] = new_cbar[i];
                });

                update_color_bar();
                update_wms(date_value);

            },start:function() {
                var color_range = 45;
                $("#cbar-slider").val(color_range);
                var date_idx = $("#slider").slider("option","value");
                var date_value = $("#select_layer option")[date_idx].value;
                cb_max = 45;

                var iter_size = cb_max / 10;
                var cbar_val = -cb_max;
                var new_cbar = [];
                for (var i=0;i<=20;i+=1){
                    new_cbar.push(parseFloat(cbar_val).toFixed(1));
                    cbar_val += iter_size;
                }
                color_bar.forEach(function(color,i){
                    color[1] = new_cbar[i];
                });

                update_color_bar();
                update_wms(date_value);
                $('#max-slider').slider("value", 45);

            }

        });

    };

    animate = function(){
        var sliderVal = $("#slider").slider("value");

        sliderInterval = setInterval(function() {
            sliderVal += 1;
            $("#slider").slider("value", sliderVal);
            if (sliderVal===layers_length - 1) sliderVal=0;
        }, animationDelay);
    };
    $(".btn-run").on("click", animate);
    //Set the slider value to the current value to start the animation at the );
    $(".btn-stop").on("click", function() {
        //Call clearInterval to stop the animation.
        clearInterval(sliderInterval);
    });

    $(".btn-increase").on("click", function() {
        clearInterval(sliderInterval);

        if(animationDelay > 250){

            animationDelay = animationDelay - 250;
            $("#speed").val((1/(animationDelay/1000)).toFixed(2));
            animate();
        }

    });

    //Decrease the slider timer when you click decrease the speed
    $(".btn-decrease").on("click", function() {
        clearInterval(sliderInterval);
        animationDelay = animationDelay + 250;
        $("#speed").val((1/(animationDelay/1000)).toFixed(2));
        animate();
    });

    add_wms = function(){
        // gs_layer_list.forEach(function(item){
        map.removeLayer(wms_layer);
        var color_str = cbar_str();
        var store_name = $("#select_layer").find('option:selected').val();
        var layer_name = 'grace:'+store_name;
        var sld_string = '<StyledLayerDescriptor version="1.0.0"><NamedLayer><Name>'+layer_name+'</Name><UserStyle><isDefault></isDefault><FeatureTypeStyle><Rule>\
        <RasterSymbolizer> \
        <ColorMap>\
        <ColorMapEntry color="#000000" quantity="'+cb_min+'" label="nodata" opacity="0.0" />'+
            color_str
            +'</ColorMap>\
        </RasterSymbolizer>\
        </Rule>\
        </FeatureTypeStyle>\
        </UserStyle>\
        </NamedLayer>\
        </StyledLayerDescriptor>';

        wms_source = new ol.source.TileWMS({
            url: wms_url,
            params: {'LAYERS':layer_name,'SLD_BODY':sld_string},
            serverType: 'geoserver',
            crossOrigin: 'Anonymous'
        });

        wms_layer = new ol.layer.Tile({
            source: wms_source
        });

        // wms_source.mergeNewParams({'SLD_BODY':sld_string});
        map.addLayer(wms_layer);

        //

    };

    update_wms = function(date_str){
        // map.removeLayer(wms_layer);
        var color_str = cbar_str();

        var layer_name = 'grace:'+date_str;
        var sld_string = '<StyledLayerDescriptor version="1.0.0"><NamedLayer><Name>'+layer_name+'</Name><UserStyle><FeatureTypeStyle><Rule>\
        <RasterSymbolizer> \
        <ColorMap> \
        <ColorMapEntry color="#000000" quantity="'+cb_min+'" label="nodata" opacity="0.0" />'+
            color_str
            +'</ColorMap>\
        </RasterSymbolizer>\
        </Rule>\
        </FeatureTypeStyle>\
        </UserStyle>\
        </NamedLayer>\
        </StyledLayerDescriptor>';

        wms_source.updateParams({'LAYERS':layer_name,'SLD_BODY':sld_string});

    };

    cbar_str = function(){
        var sld_color_string = '';
        color_bar.forEach(function(color,i){
            var color_map_entry = '<ColorMapEntry color="'+color[0]+'" quantity="'+color[1]+'" label="label'+i+'" opacity="'+color[2]+'"/>';
            sld_color_string += color_map_entry;
        });
        return sld_color_string
    };

    get_plot = function(){
        if($("#poly-lat-lon").val() == "" && $("#point-lat-lon").val() == "" && $("#shp-lat-lon").val() == ""){
            $('.warning').html('<b>No feature selected. Please create a feature using the map interaction dropdown. Plot cannot be generated without a feature.</b>');
            return false;
        }else{
            $('.warning').html('');
        }

        var datastring = $get_plot.serialize();

        $.ajax({
            type:"POST",
            url:'/apps/grace/plot-region/',
            dataType:'HTML',
            data:datastring,
            success:function(result) {
                var json_response = JSON.parse(result);
                $("#plotter").highcharts({
                    chart: {
                        type:'area',
                        zoomType: 'x',
                        height: 350
                    },
                    title: {
                        text:"Values at " +json_response.location,
                        style: {
                            fontSize: '14px'
                        }
                    },
                    xAxis: {
                        type: 'datetime',
                        labels: {
                            format: '{value:%d %b %Y}',
                            rotation: 45,
                            align: 'left'
                        },
                        title: {
                            text: 'Date'
                        }
                    },
                    yAxis: {
                        title: {
                            text: "Total Terrestrial Water Storage Anomaly (cm)"
                        }

                    },
                    exporting: {
                        enabled: true
                    },
                    series: [{
                        data:json_response.values,
                        name: "Total Terrestrial Water Storage Anomaly (cm)"
                    }]
                });
            }
        });
    };

    $("#btn-get-plot").on('click',get_plot);

    // add_test = function(){
    //
    //
    //
    //     var layer_name = 'lis:1981-01-14';
    //     var sld_string = '<StyledLayerDescriptor version="1.0.0"><NamedLayer><Name>'+layer_name+'</Name><UserStyle><FeatureTypeStyle><Rule>\
    //     <RasterSymbolizer> \
    //     <ColorMap>\
    //     <ColorMapEntry color="#000000" quantity="-0.2" label="nodata" opacity="0.0" />\
    //         <ColorMapEntry color="#1427ba" quantity="1.2" label="nodata" opacity="0.0" />\
    //         <ColorMapEntry color="#00a1ff" quantity="2.6" label="nodata" opacity="0.0" />\
    //         <ColorMapEntry color="#00ffe1" quantity="4.1" label="nodata" opacity="0.0" />\
    //         <ColorMapEntry color="#00ff48" quantity="5.5" label="nodata" opacity="0.0" />\
    //         <ColorMapEntry color="#ff0800" quantity="6.9" label="nodata" opacity="0.0" /></ColorMap>\
    //     </RasterSymbolizer>\
    //     </Rule>\
    //     </FeatureTypeStyle>\
    //     </UserStyle>\
    //     </NamedLayer>\
    //     </StyledLayerDescriptor>';
    //
    //     var lis_source = new ol.source.ImageWMS({
    //         url: 'http://127.0.0.1:8181/geoserver/wms/',
    //         params: {'LAYERS':layer_name,'SLD_BODY':sld_string},
    //         serverType: 'geoserver',
    //         crossOrigin: 'Anonymous'
    //     });
    //
    //     var lis_layer = new ol.layer.Image({
    //         source: lis_source
    //     });
    //
    //     map.addLayer(lis_layer);
    //
    // };

    /************************************************************************
     *                        DEFINE PUBLIC INTERFACE
     *************************************************************************/

    public_interface = {

    };

    /************************************************************************
     *                  INITIALIZATION / CONSTRUCTOR
     *************************************************************************/

    // Initialization: jQuery function that gets called when
    // the DOM tree finishes loading
    $(function() {
        init_vars();
        init_map();
        init_events();
        init_slider();
        gen_color_bar();

        //chart.legend.update({enabled:false});
        $maxSlider = $('#max-slider');
        $maxSlider.slider('option', 'start').call($maxSlider);

        $("#speed").val((1/(animationDelay/1000)).toFixed(2));
        $("#select_layer").change(function(){
            add_wms();
            var selected_option = $(this).find('option:selected').index();
            $("#slider").slider("value", selected_option);
        }).change();

        $("#slider").on("slidechange", function(event, ui) {
            var x = tracker[ui.value];
            chart.series[1].setData([[x,range_min],[x,range_max]]);
            var date_text = $("#select_layer option")[ui.value].text;
            $( "#grace-date" ).val(date_text); //Get the value from the slider
            var date_value = $("#select_layer option")[ui.value].value;
            update_wms(date_value);

        });
    });

    return public_interface;

}()); // End of package wrapper
// NOTE: that the call operator (open-closed parenthesis) is used to invoke the library wrapper
// function immediately after being parsed.
