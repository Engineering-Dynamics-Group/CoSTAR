%This class initially controls all input data by calling various static 
%methods of the classes DynamicalSystem, ApproxMethod and Cotinuation    
%
%The basic idea is that all errors in the input data are first analysed and
%then displayed all at the same time. This way, the user won't have to fix
%error after error but rather see them all at once. 
%
% The gatekeepers are NOT case sensitiv

classdef Gatekeeper < handle
    
    properties
        error_msg = {};                                                             %This message gets filled up by the different gate_keepers and all error messages are displayed at the end of the check.                                                            
    end

    methods
        %%%%%%%%%
        check_fields(obj,mystruct,name,needed_fieldnames,allowed_fieldnames,varargin);    	%Checks the passed mystruct if it has all needed and not more than the allowed fields. name is for display. varargin contains possible variants, if the first one fail
        
        check_data(obj,name,data,data_type,data_size,data_value);                 	        %Checks option structure entries for the correct type, the allowed values and the allowed value size.
        
        mystruct = fieldnames_to_lower(obj,mystruct);                           	        %Converts all struct fieldnames to lower case. 
        
        mystruct = fieldvalues_to_lower(obj,mystruct);                            	        %Converts all struct fields to lower case if they are a string.    

        out = my_join(obj,in);                                                    	        %gets double, a char, a cell of chars and joins them nicely to a char. (For error messages);

        speak(obj,varargin);                                                                %Displays the passed error messages

        %%%%%%%%%   
        options = m_gatekeeper(obj,options);                                                 %Main function: Does the initial security checks and calls all further gatekeepers
        
    end

end