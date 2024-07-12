% This is the help file for the tutorials
%
% The file calls the corresponding tutorial from the Tutorials folder
%
%@varargin:     optional keywords to further specify the output

function tutorial(varargin)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    default_text = [...
    '-------------------------------------------------------------\\',...
    '------------------ CoSTAR Help: Tutorials -------------------\\',...
    '-------------------------------------------------------------\\',...
    'CoSTAR provides tutorials for the user. The tutorials detailly \\' ...
    'explain specific modules of CoSTAR as well as how to use the \\' ...
    'toolbox, including the required options structure. They are \\' ...
    'recommended if you have not used CoSTAR yet. \\' ...
    ' \\' ...
    'Apart from the tutorials, CoSTAR also provides short examples, \\' ...
    'which quickly show how CoSTAR can be used. The code of the \\' ...
    'examples is identical to the code of the tutorials. However, \\' ...
    'most of the comments and explanations have been omitted. If \\' ...
    'you already know how CoSTAR works and if you just need a quick \\' ...
    'look at an example to set up your simulation, the examples are \\' ...
    'are there for you. You can go to the costarhelp page of the \\' ...
    'examples by clicking on the link below: \\']; ...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    costarhelp.s_csdisp(default_text);
    costarhelp.s_disp_fcncall('costarhelp.example()','costarhelp.example');
    fprintf('\n\n');
    fprintf('------------------------------\n');
    fprintf('\n');

    if ~isempty(varargin)
        fprintf(append('You supplied the key word(s) ',strjoin(varargin,' and '),'.\n')); 
        fprintf('\n'); 
    end

    % This function returns executable methods calls depending on the name and key words requested
    [tutorials_found,other_tutorials] = costarhelp.s_call_for_tutorial('Tutorial',costarhelp.key_words_method,costarhelp.alternative_method,costarhelp.key_words_solution_type,costarhelp.alternative_solution_type,varargin);


    switch size(tutorials_found,1)

        case 0      % This case should not exist
            
            fprintf('However and sadly, no tutorial was found. \n');  
            fprintf('Check your input or contact the support. \n');
            fprintf('\n');
            fprintf('Try to specify your call by supplying 1 or 2 key words from \n');
            fprintf('the following groups (key word group 2 is not meaningful in \n');
            fprintf('combination with EQ): \n');
            fprintf('\n');
            fprintf(append('Key word group 1: ',strjoin(costarhelp.key_words_solution_type,', '),' - or alternatively: ', strjoin(costarhelp.alternative_solution_type,', '),'.\n'));
            fprintf(append('Key word group 2: ',strjoin(costarhelp.key_words_method,', '),' - or alternatively: ', strjoin(costarhelp.alternative_method,', '),'.\n'));
            fprintf('\n');
            fprintf('Example: --> costarhelp.tutorial("PS","FGM") \n');
            fprintf('\n');
            fprintf('-------------------------------------------------------------\n');
            fprintf('\n\n\n');


        case 1      % One specific file was found. The specified file is then opened
            
            edit(tutorials_found{1})

            fprintf('Your requested tutorial has been opened in the Matlab editor.\n');

            if ~isempty(other_tutorials)
                fprintf('\n');
                fprintf('------------------------------\n');
                fprintf('\n');
                costarhelp.s_csdisp(['There are other, related tutorials, which you might find \\' ...
                    'interesting. You can open them in the Matlab editor by \\clicking on the links below:\\ '])

                for k = 1:size(other_tutorials,1)
                    fprintf(append('<a href="matlab: edit(''',other_tutorials{k},''')">',other_tutorials{k},'</a>\n'));
                end
            end

            fprintf('\n');
            fprintf('-------------------------------------------------------------\n');
            fprintf('\n\n\n');


        otherwise   % More than one specific file was found. Hyperlinks to the found files are displayed
            
            fprintf(['There are multiple tutorials for your request. You can open \n' ...
                     'them in the Matlab editor by clicking on the links below:\n']); 
            fprintf('\n');
            for k = 1:size(tutorials_found,1)
                fprintf(append('<a href="matlab: edit(''',tutorials_found{k},''')">',tutorials_found{k},'</a>\n'));
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
            fprintf('Example: --> costarhelp.tutorial("PS","FGM") \n');
            fprintf('\n');
            fprintf('-------------------------------------------------------------\n');
            fprintf('\n\n\n');

    end
    
end