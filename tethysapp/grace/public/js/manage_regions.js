/*****************************************************************************
 * FILE:    GRACE MANAGE REGIONS
 * DATE:    19 MAY 2017
 * AUTHOR: Sarva Pulla
 * COPYRIGHT: (c) Brigham Young University 2017
 * LICENSE: BSD 2-Clause
 *****************************************************************************/


/*****************************************************************************
 *                      LIBRARY WRAPPER
 *****************************************************************************/

var GRACE_MANAGE_REGIONS = (function() {
    // Wrap the library in a package function
    "use strict"; // And enable strict mode for this library

    /************************************************************************
    *                      MODULE LEVEL / GLOBAL VARIABLES
    *************************************************************************/
    var m_uploading_data, m_results_per_page;

    /************************************************************************
     *                    PRIVATE FUNCTION DECLARATIONS
     *************************************************************************/
    var initializeTableFunctions, getTablePage, displayResultsText;


    /************************************************************************
     *                    PRIVATE FUNCTION IMPLEMENTATIONS
     *************************************************************************/

    m_results_per_page = 5;

    initializeTableFunctions = function() {

        //handle the submit delete event
        $('.submit-delete-region').click(function(){
        var data = {
            region_id: $(this).parent().parent().parent()
                                .find('.region-id').text()
                };
            //update database
            var xhr = deleteRowData($(this), data);
            if (xhr != null) {
                xhr.done(function (data) {
                    if ('success' in data) {
                        addSuccessMessage("Region Successfully Deleted!");
                        var num_regions_data = $('#manage_regions_table').data('num_regions');
                        var page = parseInt($('#manage_regions_table').data('page'));
                        $('#manage_regions_table').data('num_regions', Math.max(0, parseInt(num_regions_data) - 1));
                        if (parseInt($('#manage_regions_table').data('num_regions')) <= m_results_per_page * page) {
                            $('#manage_regions_table').data('page', Math.max(0, page - 1));
                        }
                        getTablePage();
                    }
                });
            }
        });

        displayResultsText();
        if (m_results_per_page >= $('#manage_regions_table').data('num_regions')) {
            $('[name="prev_button"]').addClass('hidden');
            $('[name="next_button"]').addClass('hidden');
        }

        //pagination next and previous button update
        $('[name="prev_button"]').click(function(){
            var page = parseInt($('#manage_regions_table').data('page'));
            $('#manage_regions_table').data('page', Math.max(0, page-1));
            getTablePage();
        });
        $('[name="next_button"]').click(function(){
            var page = parseInt($('#manage_regions_table').data('page'));
            $('#manage_regions_table').data('page', Math.min(page+1,
                                                Math.floor(parseInt($('#manage_regions_table').data('num_regions')) / m_results_per_page - 0.1)));
            getTablePage();
        });
    };


    displayResultsText = function() {
        //dynamically show table results display info text on page
        var page = parseInt($('#manage_regions_table').data('page'));
        var num_regions_data = $('#manage_regions_table').data('num_regions');
        var display_min;
        if (num_regions_data == 0){
            display_min = 0
        }
        else{
            display_min = ((page + 1) * m_results_per_page) - (m_results_per_page - 1);
        }
        var display_max = Math.min(num_regions_data, ((page + 1) * m_results_per_page));
        $('[name="prev_button"]').removeClass('hidden');
        $('[name="next_button"]').removeClass('hidden');
        if (page == 0){
            $('[name="prev_button"]').addClass('hidden');
        } else if (page == Math.floor(num_regions_data / m_results_per_page - 0.1)) {
            $('[name="next_button"]').addClass('hidden');
        }
        if (num_regions_data != 0) {
            $('#display-info').append('Displaying regions ' + display_min + ' - ' +
                display_max + ' of ' + num_regions_data);
        }else {
            $('#display-info').append('No regions to display' + '<br>To add one, ' +
                'click <a href="../add-region">here</a>.');
        }
    };

    getTablePage = function() {
        $.ajax({
            url: 'table',
            method: 'GET',
            data: {'page': $('#manage_regions_table').data('page')},
            success: function(data) {
                $("#manage_regions_table").html(data);
                initializeTableFunctions();
            }
        });
    };


    /************************************************************************
    *                  INITIALIZATION / CONSTRUCTOR
    *************************************************************************/

    $(function() {
        m_uploading_data = false;
        getTablePage();
    }); //document ready
}()); // End of package wrapper