/*****************************************************************************
 * FILE:    GRACE ADD GEOSERVER
 * DATE:    18 MAY 2017
 * AUTHOR: Sarva Pulla
 * COPYRIGHT: (c) Brigham Young University 2017
 * LICENSE: BSD 2-Clause
 *****************************************************************************/

/*****************************************************************************
 *                      LIBRARY WRAPPER
 *****************************************************************************/

var GRACE_ADD_GEOSERVER = (function() {
    // Wrap the library in a package function
    "use strict"; // And enable strict mode for this library

    /************************************************************************
     *                      MODULE LEVEL / GLOBAL VARIABLES
     *************************************************************************/
    var public_interface;				// Object returned by the module



    /************************************************************************
     *                    PRIVATE FUNCTION DECLARATIONS
     *************************************************************************/

    var add_geoserver,init_jquery,reset_alert,reset_form;

    /************************************************************************
     *                    PRIVATE FUNCTION IMPLEMENTATIONS
     *************************************************************************/

    init_jquery = function(){
    };

    //Reset the alerts if everything is going well
    reset_alert = function(){
        $("#message").addClass('hidden');
        $("#message").empty()
            .addClass('hidden')
            .removeClass('alert-success')
            .removeClass('alert-info')
            .removeClass('alert-warning')
            .removeClass('alert-danger');
    };

    //Reset the form when the request is made succesfully
    reset_form = function(result){
        if("success" in result){
            $("#geoserver-name-input").val('');
            $("#geoserver-url-input").val('');
            $("#geoserver-username-input").val('');
            $("#geoserver-password-input").val('');
            addSuccessMessage('Geoserver Upload Complete!');
        }
    };

    add_geoserver = function(){
        reset_alert();
        var geoserver_name = $("#geoserver-name-input").val();
        var geoserver_url = $("#geoserver-url-input").val();

        var geoserver_username = $("#geoserver-username-input").val();
        var geoserver_password = $("#geoserver-password-input").val();

        if(geoserver_name == ""){
            addErrorMessage("Geoserver Name cannot be empty!");
            return false;
        }else{
            reset_alert();
        }
        if(geoserver_url == ""){
            addErrorMessage("Geoserver Url cannot be empty!");
            return false;
        }else{
            reset_alert();
        }
        if(geoserver_username == ""){
            addErrorMessage("Geoserver Username cannot be empty!");
            return false;
        }else{
            reset_alert();
        }
        if(geoserver_password == ""){
            addErrorMessage("Geoserver Password cannot be empty!");
            return false;
        }else{
            reset_alert();
        }

        if(geoserver_url.includes('/geoserver/rest') == false){
            addErrorMessage("Geoserver should have /geoserver/rest at the end of the url");
            return false;
        }else{
            reset_alert();
        }

        if (geoserver_url.substr(-1) !== "/") {
            geoserver_url = geoserver_url.concat("/");
        }
        var data = {"geoserver_name":geoserver_name,"geoserver_url":geoserver_url,"geoserver_username":geoserver_username,"geoserver_password":geoserver_password};

        var xhr = ajax_update_database("submit",data);
        xhr.done(function(return_data){
            if("success" in return_data){
                reset_form(return_data);
            }else if("error" in return_data){
                addErrorMessage(return_data["error"]);
            }
        });



    };
    $("#submit-add-geoserver").click(add_geoserver);

    /************************************************************************
     *                        DEFINE PUBLIC INTERFACE
     *************************************************************************/
    /*
     * Library object that contains public facing functions of the package.
     * This is the object that is returned by the library wrapper function.
     * See below.
     * NOTE: The functions in the public interface have access to the private
     * functions of the library because of JavaScript function scope.
     */
    public_interface = {

    };

    /************************************************************************
     *                  INITIALIZATION / CONSTRUCTOR
     *************************************************************************/

    // Initialization: jQuery function that gets called when
    // the DOM tree finishes loading
    $(function() {
        init_jquery();

    });

    return public_interface;

}()); // End of package wrapper
// NOTE: that the call operator (open-closed parenthesis) is used to invoke the library wrapper
// function immediately after being parsed.