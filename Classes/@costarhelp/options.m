% This is the help method for the option system.
% It displays the help message on the command window and gives further options

function options(varargin)
        

    fprintf('\n');
    if ~isempty(varargin)
        disp('Attention: costarhelp.options allows no further specifying arguments!')
    end
    fprintf('\n');

    costarhelp.s_csdisp([...
    '-------------------------------------------------------------\\',...
    '------------- CoSTAR Help for options structure -------------\\',...
    '-------------------------------------------------------------\\',...
    'In order to compute a single solution point or to continue \\' ...
    'a solution, you can call CoSTAR by: \\' ...
    ' \\' ...
    ' --> [S,DYN] = costar(options);\\']); ...
    fprintf('\n');
    disp('Here, the structure options contains all information needed for ');
    disp('the simulation. options itself consists of different subordinated ');
    disp('structures. You can create these structures like a normal struct ');
    disp('in Matlab, but you have to use the costaropts() function. E.g. ');
    disp('for the options.system structure:');
    fprintf('\n');
    fprintf('--> options.system = costaropts("order",1,"rhs",FCNhandle,"dim",2) \n');
    fprintf('\n');
    fprintf('Click the links below to learn more!\n');
    
    fprintf('\n');
    disp('Mandatory costaropts structs:');
    fprintf('\n');

    costarhelp.s_disp_fcncall('costarhelp.system()','options.system');     fprintf(':            defines the system to be continued \n');
    costarhelp.s_disp_fcncall('costarhelp.opt_sol()','options.opt_sol');   fprintf([':           contains info on the solution type \n',...
                                                                     '                           branch and method \n']);
    costarhelp.s_disp_fcncall('costarhelp.opt_init()','options.opt_init'); fprintf([':          defines the initial condition for searching \n',...
                                                                    '                           the first branch point (approx method specific) \n']);
    
    fprintf('\n');
    costarhelp.s_csdisp('Optional or mandatory costaropts structs:') 
    fprintf('\n');
    
    costarhelp.s_disp_fcncall('costarhelp.opt_cont()','options.opt_cont');                     fprintf(':            defines information for the branch continuation \n');
    costarhelp.s_disp_fcncall('costarhelp.opt_approx_method()','options.opt_approx_method');  fprintf([':   sepcifies information for the chosen \n',...
                                                                              '                             approximation method (exception: equilibrium) \n']);
    costarhelp.s_disp_fcncall('costarhelp.opt_stability()','options.opt_stability');            fprintf([':       defines infos for stability and bifuractions \n',...
                                                                                    '                             (solution type specific) \n']);

    fprintf('\n\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('\n\n\n');

end