%This is the help file for the opt_init structure
%
%The file calls the corresponding s_help_xxx static methods from the corresponding classes
%
%@varargin:     optional keywords to further specify the output

function opt_init(varargin)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    default_text = [...
    '-------------------------------------------------------------\n',...
    '----------------- CoSTAR Help for opt_init  -----------------\n',...
    '-------------------------------------------------------------\n',...
    'The opt_init costaropts struct is used to define initial\n',...
    'conditions to compute the first curve point / solution for\n' ...
    'one of the following approximation methods:       \n',...
     strjoin(costarhelp.alternative_method,' or '),'. \n' ...
     '\n' ...
     '------------------------------\n' ...
     '\n'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf(default_text);

    % Call for help! This function returns executable methods calls for the help methods in the different class folders corresponding to the call.
    [help_method_name,base_method_names] = costarhelp.s_call_for_help('opt_init',costarhelp.key_words_method,costarhelp.alternative_method,costarhelp.key_words_solution_type,costarhelp.alternative_solution_type,varargin);
    

    switch size(help_method_name,1)
    
        case 0      % This case should not exist
            
            if ~isempty(varargin)
                fprintf(append('You supplied the key word(s) ',strjoin(varargin,' and '),'.\n')); 
                fprintf('\n');
            end
            fprintf('However and sadly, no help file was found. \n');  
            fprintf('Check your input or contact the support. \n');
            fprintf('\n');
            fprintf('Try to specify your call by supplying 1 or 2 key words from \n');
            fprintf('the following groups (key word group 2 is not meaningful in \n');
            fprintf('combination with EQ): \n');
            fprintf('\n');
            fprintf(append('Key word group 1: ',strjoin(costarhelp.key_words_solution_type,', '),' - or alternatively: ', strjoin(costarhelp.alternative_solution_type,', '),'.\n'));
            fprintf(append('Key word group 2: ',strjoin(costarhelp.key_words_method,', '),' - or alternatively: ', strjoin(costarhelp.alternative_method,', '),'.\n'));
            fprintf('\n');
            fprintf('Example: --> costarhelp.opt_init("PS","fourier-galerkin") \n');
            fprintf('\n');
            fprintf('-------------------------------------------------------------\n');
            fprintf('\n\n\n');


        case 1      % One specific help file was found. The specified help file is then called.
            
            if ~isempty(varargin)
                fprintf(append('You supplied the key word(s) ',strjoin(varargin,' and '),'.\n')); 
            end

            help_struct = eval(help_method_name{1,1});          %Get the help struct.
            costarhelp.s_disp_help_text(help_struct);           %Display the help structure

            if ~isempty(base_method_names)
                fprintf('\n');
                costarhelp.s_csdisp(['There are other, related help files, which you might find \\' ...
                    'interesting. Choose one of the following options by \\clicking on the link:\\ '])

                for k = 1:size(base_method_names,1)
                    fprintf(append('<a href="matlab:','costarhelp.s_disp_help_text(',base_method_names{k,1},')','">','Help for',' ',base_method_names{k,2},'</a>'));
                    fprintf('\n');
                end
            end

            fprintf('\n\n');
            fprintf('-------------------------------------------------------------\n');
            fprintf('\n\n\n');


        otherwise   % More than one specific help file was found. Hyperlinks to the found help files are displayed.

            if ~isempty(varargin)
                fprintf(append('You supplied the key word(s) ',strjoin(varargin,' and '),'.\n')); 
                fprintf('\n');
            end
            fprintf(['There are multiple help files for your help request for \n' ...
                     'opt_init. Choose one of the following options by clicking \n' ...
                     'on the link:\n']);     
            fprintf('\n');
            for k = 1:size(help_method_name,1)
                fprintf(append('<a href="matlab:','costarhelp.s_disp_help_text(',help_method_name{k,1},')','">','Help for',' ',help_method_name{k,2},'</a>'));
                fprintf('\n');
            end

            fprintf('\n');
            fprintf('------------------------------\n');
            fprintf('\n');
            fprintf('Alternatively, you can specify your call by supplying 1 or 2 \n')
            fprintf('key words from the following groups (key word group 2 is not \n')
            fprintf('meaningful in combination with EQ): \n');
            fprintf('\n');
            fprintf(append('Key word group 1: ',strjoin(costarhelp.key_words_solution_type,', '),' - or alternatively: ', strjoin(costarhelp.alternative_solution_type,', '),'.\n'));
            fprintf(append('Key word group 2: ',strjoin(costarhelp.key_words_method,', '),' - or alternatively: ', strjoin(costarhelp.alternative_method,', '),'.\n'));
            fprintf('\n');
            fprintf('Example: --> costarhelp.opt_init("PS","fourier-galerkin") \n');
            fprintf('\n');
            fprintf('-------------------------------------------------------------\n');
            fprintf('\n\n\n');
    
    end

end