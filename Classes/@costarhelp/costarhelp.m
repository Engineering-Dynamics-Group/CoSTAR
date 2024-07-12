%This is a class for better organizing the help files in CoSTAR.
%
%@varargin:
%   empty:      general information on costar and the costarhelp tool
% 1 entry:      either an option name or a function name (this will call the help file of the corresponding function)
% 2 entries:    one must be an option or postprocessing command, the second can sepcify the information, e.g. EQ, PS, QPS or FGM, Shooting, etc.
% 3 entries:    one must be an option or postprocessing command, the second and third can sepcify the information further, e.g. EQ, PS, QPS or FGM, Shooting, etc.



classdef costarhelp < handle

    
    properties(Constant)

            %% These are the key words for focussing the help request
            key_words_method          = {'FDM','FGM','SHM'};                                    %Possible words for the first or second input; Second half are alternative words
            alternative_method        = {'finite-difference','fourier-galerkin','shooting'};    %Alternative words for the possible_key_words_method: ATTENTION: Same order!
            key_words_solution_type   = {'EQ','PS','QPS'};                                      %Possible words for the other input; Second half are alternative words
            alternative_solution_type = {'equilibrium','periodic','quasi-periodic'};            %Alternative words for the possible_key_words_solution_type: ATTENTION: Same order!
            key_words_stability = {'SHM'};
            alternative_words_stability = {'shooting'};



    end


    methods
        %% Constructor

        function obj = costarhelp()  

                costarhelp.costar();       %Display the main help function for costar.

        end


    end

    methods(Static)
        
        %% Process Methods
        [help_method_name,base_method_names] = s_call_for_help(help_name,possible_words_1,alternative_1,possible_words_2,alternative_2,varargin);   %Checks the input data and searches for the corresponding help files in the Class folder
        [files_found,other_files] = s_call_for_tutorial(name,possible_words_1,alternative_1,possible_words_2,alternative_2,varargin);      %Checks the input data and searches for the corresponding tutorial or minimal example files
        s_disp_help_text(s_help_method);                                                                %Display either a help text or the help option struct
        s_csdisp(mystring);                                                                             %This displays the string given to the function WHILE the Latex line break command "\\" can be used for line breaks
        s_disp_fcncall(fcn,name);                                                                       %This function displays a clickable hyperlink called name, which executes the function fcn
        s_disp_opt_help_struct(help_struct);                                                            %This function displays the option help structure, defined in the corresponding s_help_opt_xxx files
              
        formatted_str = s_format_string(str, max_width);                                                %This function formats a string to a given width and inserts line breaks
        output = s_format_string_to_cell(str, width,breaktype);                                         %This function formats a string to a given width and and returns cell arrayss

        %% Help functions: These purely display help information on the command window
        costar(varargin);                   %Main help message for costar
        options(varargin);                  %Overview help message for the different option structures.

        system(varargin);                   %Displays the help message for the system option structure
        opt_init(varargin);                 %Displays the help messages for the opt_init option structure  
        opt_sol(varargin);                  %Displays the help messages for the opt_sol option structure  

        opt_cont(varargin);                 %Displays the help messages for the opt_cont option structure     
        opt_approx_method(varargin);        %Displays the help messages for the opt_approx_method option structure   
        opt_stability(varargin);            %Displays the help messages for the opt_stability option structure 

        solget(varargin);                   %Displays the help messages for the solget option structure  
        solplot(varargin);                  %Displays the help messages for the solplot option structure  
        contplot(varargin);                 %Displays the help messages for the contplot option structure  
    
        example(varargin);                  %Displays info text for the minimal examples
        tutorial(varargin);                 %Displays info text for the tutorials
    end






end