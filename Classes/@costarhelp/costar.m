%This is the main help file for the costarhelp function
%
%@obj:      CShelp class object

function costar()
    
    disp('-------------------------------------------------------------');
    disp('------------------------ CoSTAR Help ------------------------');
    disp('-------------------------------------------------------------');
    fprintf('\n');
    fprintf(['CoSTAR is a Matlab path continuation toolbox written by the \n', ...
        '<a href="https://www.uni-kassel.de/maschinenbau/en/institute/mechanik/fachgebiete/engineering-dynamics/frontpage">Engineering Dynamics Group of the University of Kassel, Germany.</a> \n ']);
    fprintf('\n');
    fprintf(['CoSTAR enables you to compute single stationary solutions or \n', ...
             'to continue a  branch of equilibria, periodic or quasi-periodic \n', ...
             'solutions with different numerical techniques, whilst tracking \n',...
             'the solution stability and identify occurring bifurcations. \n']);
    

    fprintf('\n\n');

    disp('-------------------------------------------------------------');
    disp('------ Single Solution points or Branch Continuation  -------');
    disp('-------------------------------------------------------------');
    fprintf('\n');
    fprintf(['In order to compute a single solution point or to continue \n',...
             'a solution, you can call CoSTAR by \n' ...
              '\n'...
               ' --> [S,DYN] = costar(options); \n' ...
               '\n' ...
               'This call is always identical. \n']);
    disp('Here, the structure options contains all information needed for ');
    disp('the simulation. options itself consists of different subordinated ');
    disp('structures. You can create these structures like a normal struct ');
    disp('in Matlab, but you have to use the costaropts() function. E.g. ');
    disp('for the options.system structure:');
    fprintf('\n');
    fprintf('--> options.system = costaropts("order",1,"rhs",FCNhandle,"dim",2) \n');

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
    costarhelp.s_disp_fcncall('costarhelp.opt_approx_method()','options.opt_approx_method');   fprintf([':   sepcifies information for the chosen \n',...
                                                                              '                             approximation method (exception: equilibrium) \n']);
    costarhelp.s_disp_fcncall('costarhelp.opt_stability()','options.opt_stability');            fprintf([':       defines infos for stability and bifuractions \n',...
                                                                                    '                             (solution type specific) \n']);


    fprintf('\n\n');

    disp('-------------------------------------------------------------');
    disp('---------------------- Post-Processing ----------------------');
    disp('-------------------------------------------------------------');
    fprintf('\n');
    fprintf('The Solution object S, received by a successful costar call, \n');
    fprintf('offers multiple postprocessing methods. You can call one of \n');
    fprintf('the following three methods: \n');
    fprintf('\n');

    costarhelp.s_disp_fcncall('costarhelp.solget()'  ,'S.solget(DYN,solget_opts)');    fprintf(':       extracts solution data defined in solget_opts \n');
    costarhelp.s_disp_fcncall('costarhelp.solplot()' ,'S.solplot(DYN,solplot_opts)');    fprintf([':     displays the solution plot at a single \n',...
                                                                     '                                 solution point, defined in solplot_opts \n']);
    costarhelp.s_disp_fcncall('costarhelp.contplot()','S.contplot(DYN,contplot_opts)'); fprintf([':   displays a continuation plot of the \n',...
                                                                    '                                 entire curve defined in contplot_opts \n']);

    fprintf('\n');
    disp('solget_opts, solplot_opts and contplot_opts are all costaropts objects.');


    fprintf('\n\n');
    disp('-------------------------------------------------------------');
    disp('----------------------- Further Help  -----------------------');
    disp('-------------------------------------------------------------');  
    fprintf('\n');
    fprintf('In order to find out more on the costaropts objects, call the \ncostarhelp function:\n \n');
    costarhelp.s_csdisp(['Example: --> costarhelp.opt_approx_method \\' ...
        ' \\' ...
        'This displays all help massages available for opt_approx_method \\' ...
        'for all solutions types and methods. \\ \\']);
    costarhelp.s_csdisp('Help files are available for:'); ...

    costarhelp.s_disp_fcncall('costarhelp.costar()'             ,'costarhelp.costar');                              fprintf('\n');
    costarhelp.s_disp_fcncall('costarhelp.options()'            ,'costarhelp.options');                             fprintf('\n \n');

    costarhelp.s_disp_fcncall('costarhelp.system()'             ,'costarhelp.system');                              fprintf('\n');
    costarhelp.s_disp_fcncall('costarhelp.opt_sol()'            ,'costarhelp.opt_sol');                             fprintf('\n');
    costarhelp.s_disp_fcncall('costarhelp.opt_init()'           ,'costarhelp.opt_init("Key_1","Key_2")');           fprintf('\n');
    costarhelp.s_disp_fcncall('costarhelp.opt_cont()'           ,'costarhelp.opt_cont');                            fprintf('\n');
    costarhelp.s_disp_fcncall('costarhelp.opt_approx_method()'  ,'costarhelp.opt_approx_method("Key_1","Key_2")');  fprintf('\n');
    costarhelp.s_disp_fcncall('costarhelp.opt_stability()'      ,'costarhelp.opt_stability("Key_1","Key_2")');      fprintf('\n \n');

    costarhelp.s_disp_fcncall('costarhelp.solget()'             ,'costarhelp.solget');                              fprintf('\n');
    costarhelp.s_disp_fcncall('costarhelp.solplot()'            ,'costarhelp.solplot');                             fprintf('\n');
    costarhelp.s_disp_fcncall('costarhelp.contplot()'           ,'costarhelp.contplot');                            fprintf('\n \n');

    costarhelp.s_disp_fcncall('costarhelp.example()'            ,'costarhelp.example("Key_1","Key_2")');            fprintf('\n');
    costarhelp.s_disp_fcncall('costarhelp.tutorial()'           ,'costarhelp.tutorial("Key_1","Key_2")');           fprintf('\n \n');
    

    costarhelp.s_csdisp(['In order to make your call more specific, you can add up to 2 keys:\\' ...
        ' \\' ...
        'Example: --> costarhelp.opt_approx_method("fourier-galerkin","periodic") \\' ...
        ' \\' ...
        'This displays the help file for the opt_approx_method for the \\Fourier-Galerkin Method for periodic solutions.']);

    fprintf('You can choose from the following (not every combination available):\n');
    fprintf('\n');
    fprintf(append('Key word group 1: ',strjoin(costarhelp.key_words_solution_type,', '),' - or alternatively ', strjoin(costarhelp.alternative_solution_type,', '),'.\n'));
    fprintf(append('Key word group 2: ',strjoin(costarhelp.key_words_method,', '),' - or alternatively ', strjoin(costarhelp.alternative_method,', '),'.\n'));
    fprintf('\n');
    costarhelp.s_csdisp(['For further information contact the staff of the \\'...
             '<a href="https://www.uni-kassel.de/go/technische-dynamik">Engineering Dynamics Group of the University of Kassel, Germany.</a>']);
    fprintf('\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('\n\n\n');
   
end