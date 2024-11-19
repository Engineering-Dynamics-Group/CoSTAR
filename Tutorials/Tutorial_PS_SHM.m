%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      CoSTAR                       %
%   Continuation of Solution Torus AppRoximations   %
%                                                   %
%                    Tutorial:                      %
%               Periodic Solutions                  %
%         - (Multiple) Shooting Method -            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%           Welcome to the CoSTAR Tutorials!        %
%
% The tutorials comprehensively explain certain CoSTAR modules, which is why they are the perfect starting point for CoSTAR beginners.
% They are highly recommended if you are not yet familiar with the CoSTAR toolbox. 
% 
% If you have already used CoSTAR, you may find the CoSTAR examples helpful.
% These provide example code to briefly show how the toolbox can be set up and how a certain CoSTAR module can be used.
% There is a corresponding example for each tutorial and the code of both of them is identical.
%
% In this tutorial, we will look at how to compute periodic solutions using the (Multiple) Shooting Method (SHM) to approximate the ...           
% solution (associated example: Example_PS_SHM). We will learn how to set up CoSTAR and how to correctly define all the required ...
% settings on the basis of two different examples:
%
%  1. Continuation of the Duffing oscillator with external forcing (non-autonomous system)
%       1.1 Introduction
%       1.2 Continuation of excitation frequency
%           1.2.1 Setting useful variables
%           1.2.2 The CoSTAR settings
%           1.2.3 Calling CoSTAR and running the simulation
%       1.3 Postprocessing: Plot of time-dependent behaviour of solution x(tau,mu)
%       1.4 Change of continuation parameter from excitation frequency to coefficient of nonlinear stiffness
%       ->  Usage of all SHM specific fields of "options.opt_init" and "options.opt_approx_method"
%
%  2. Continuation of the van der Pol oscillator (autonomous system)
%       2.1 Introduction
%       2.2 Continuation of nonlinear damping
%           2.2.1 Setting useful variables
%           2.2.2 The CoSTAR settings
%           2.2.3 Calling CoSTAR and running the simulation
%       2.3 Postprocessing: Plot of time-dependent behaviour of solution x(tau,mu)
%       ->  Usage of all SHM specific fields of "options.opt_init" and "options.opt_approx_method"
%
% NOTE: Both examples can be worked through independently.
%
% IMPORTANT: Explanations can be found at the beginning of new lines as well as right next to code, so make sure that you do not ...
%            miss some helpful comments!
%            Furthermore, it is advised to execute the code line by line or in blocks when going through this tutorial and do not ...
%            execute the complete script by clicking "Run" or pressing "F5". To do that, you can select the desired lines of code ...
%            to be executed and press "F9" or you can click right on the selected code and choose "Evaluate Selection in Command ...
%            Window". This should help to better understand what is going on and how CoSTAR works. 
%
% Alright, let's clean our "desk" and let's get started!
clear variables; clc; close all;                                        % clear workspace; clear command window; close all figures

% We need to add the main CoSTAR folder (and all subfolders) to MATLAB's search path. 
% Assuming you are running this script within the "Tutorials" subfolder of CoSTAR, this is done by:
% addpath(genpath('..\'))                                               % genpath() generates a search path containing all subfolders
% (The command is set as comment in order to prevent unwanted folders to be added to MATLAB's search path.) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Example No. 1: Duffing Oscillator        %%
%               Non-Autonomous System               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%               1.1  Introduction                   %
%
% The first example we have a look on features the Duffing equation with external forcing. The dimensionless form of the ...
% equation reads:     x'' + 2*D*x' + x + kappa*x^3 = g*cos(eta*tau)     (1)
% where x is the deflection, D is the damping factor, kappa is the coefficient of the non-linear stiffness, g is the amplitude ...
% of the excitation, eta is the dimensionless excitation frequency and tau is the dimensionless time. Following that, ...
% x' = dx/dtau and x'' = dx'/dtau holds.
% However, CoSTAR always needs the equation to be of the form dz/dt = f(t,z,param) (param represents all additional parameters), ...
% which is why we have to transform (1) into a system of two differential equations of first order. To do that, we introduce the ...
% state variable vector z = [z1; z2] = [x; x'], leading to:
% dz/dtau = f(tau,z,kappa,D,eta,g)  (2)  <=>  z1' = z2
%                                             z2' = -2*D*z2 - z1 - kappa*z1^3 + g*cos(eta*tau) 
%
% Our aim is to approximate the solution of equation (2) z(tau) for eta in the range of [0.01, 2.5] using the SHM. This problem ...
% can be rewritten to finding the zero set of the function g: R^(n+1) -> R^n with g = 0 and for eta in [0.01, 2.5], which ...
% defines a curve in (n+1)-dimensional space (n equals the dimension of the system. Here, n = 2 holds). 
% Hence, a solution z(tau) of (2) is also denoted as a "point on the curve" in the following explanations. CoSTAR will approximate ...
% discrete points on this curve by executing a path continuation. As the frequency eta is varied, the continuation parameter is ...
% mu = eta. The remaining parameters in equation (2) are set to:
D = 0.05;
kappa = 0.3;
g = 1;


%     1.2 Continuation of excitation frequency      %
%
%          1.2.1 Setting useful variables           
%
% Before we specify the settings needed by CoSTAR, we introduce some important variables, which we will use later on. This is not ...
% necessarily needed but it helps to keep an overview of the most important variables/settings. Furthermore, important aspects of ...
% CoSTAR are explained here, making the following section "1.2.2 The CoSTAR settings" clearer.
mu_limit = [0.01, 2.5];                 % mu_limits = eta_limits = [0.01, 2.5] sets the upper and bottom limit of our continuation.
%
eta0 = mu_limit(1);                     % eta0 defines the eta-value at which the continuation begins.
%
Fcn =  @(t,z,param) duffing_ap(t,z,param);  % Fcn is a function handle containing the right-hand side (RHS) of equation (2). 
% You may have noticed that Fcn is a function of t, z and param, whereas the right-hand side of (2) is a function of tau, z, ...
% kappa, D, eta and g. This is because CoSTAR ALWAYS requires the definition of Fcn as "Fcn = @(t,z,param) ..." when calculating ...
% periodic solutions. param is a cell array storing all parameters, apart from t and z, that need to be passed to (2). The time t ...
% as well as the state variable vector z are always treated separately and therefore not included in param. Here, we need to define ...
% param since we need to pass the parameters kappa, D, eta and g to Fcn. The name of the "time argument" being t (and not tau) is ...
% required by CoSTAR. It is equal to our dimensionless time tau.
%
% TASK: Open the function "duffing_ap" (click the right mouse button right on it and choose 'Open "duffing_ap"' or go into the "RHS" ...
%       subfolder located in the main path of CoSTAR and double click on the function). Have a look at how the equation is defined ...
%       and compare that to our system right next to equation (2).
% IMPORTANT: Please pay attention to how the state variables z1 and z2 are defined and how the parameters D, kappa, eta and g ...
%            are extracted from param.
%
% Now we define the "param" array. It is important that it is a cell array AND the order of the parameters placed in param must ...
% correspond to the definitions in the "duffing_ap" function! For example, kappa is defined by kappa = param{2} in the "duffing_ap" ...
% function. Therefore, kappa has to be the SECOND variable in the "param" array. For the variable eta, we use eta0 since the ...
% continuation starts at the eta-value of eta0.
param = {kappa, D, eta0, g};   
% Next, we have to tell CoSTAR where the continuation parameter is located within param. We also call this the "active parameter". ...
% This is done by defining an integer variable which is used to get the continuation parameter from param by cell indexing, i.e. ...
% continuation parameter mu = param{active_parameter}. As we want to continue the curve via eta and the corresponding value eta0 ...
% is the THIRD element of param, we set:
active_parameter = 3;
%
% Going on, we need to set the non-autonomous frequency of the external forcing (eta). This sounds a bit confusing since eta is the ...
% continuation parameter mu and we already set the range of the continuation by defining mu_limit. So why do we need to specify the ...
% non-autonomous frequency again? The reason is that eta CAN be the continuation parameter but does NOT HAVE TO BE the continuation ...
% parameter. In order to deal with this uncertainty, we introduce the variable "non_auto_freq" which sets the non-autonomous ...
% frequency. In this case, the non-autonomous frequency is the continuation parameter mu. Due to this dependency, we have to set ...
% "non_auto_freq" as function handle @(mu):
non_auto_freq = @(mu) mu;
% "non_auto_freq" ALWAYS has to be a function handle @(mu) if it depends on mu. This also allows a great flexibility since the ...
% non-autonomous frequency can be an arbitrary function of mu, e.g. @(mu) 10*mu. In the other case that the non-autonomous ...
% frequency is independent from mu and therefore constant throughout the continuation (i.e. another parameter is the continuation ...
% parameter), we can define "non_auto_freq" either as double (e.g. non_auto_freq = 1.5) or as function handle ...
% (e.g. non_auto_freq = @(mu) 1.5).
%
% Last but not least, we set a variable called "initial condition" (IC). Using the SHM, the initial condition defines a point in ...
% state space where the nonlinear system solver fsolve starts to find the first approximate solution of equation (2) at mu0.
% Remember that the SHM iterates towards a point z_po on the periodic orbit (which is a closed trajectory) of sought solution ...
% z(tau,mu) in state space. The time-dependent approximate solution z(tau,mu) is obtained by time integration starting at z_po.
% The initial condition should be as close as possible to the periodic orbit.
IC = [g; 0];


%            1.2.2 The CoSTAR settings              
%
% All information and settings which CoSTAR needs are defined in a structure array called "options". 
% The structure "options" itself consists of structure arrays that define important parameters used by different modules of CoSTAR.
% 
% CoSTAR ALWAYS needs the following "option" structures: "options.system", "options.opt_sol" and "options.opt_init".
% Depending on the system and problem, CoSTAR might need some or all of the following "option" structures additionally: 
% "options.opt_approx_method":  Required when calculating (quasi-) periodic solutions.
% "options.opt_cont":           Required when running a continuation.
% "options.opt_stability":      Optional when calculating the stability of solutions.
% In this case, "options.opt_approx_method" and "options.opt_cont" are additionally needed.
% 
% Each of the "option" structures contain mandatory fields which ALWAYS have to be set. 
% Furthermore, there are optional fields that CAN be set. Most of the optional fields have a default value. 
% However, even optional fields NEED to be set in some cases (denoted as "Opt-Need.", e.g. 'param' in this example, see below).
% 
% IMPORTANT: All option structures must be created by the "costaropts" function. The syntax of the arguments of "costaropts" is ...
%            equivalent to the MATLAB struct function, i.e. costaropts(fieldname1,value1,...,fieldnameN,valueN).
%
% Let us start with the "options.system" structure. There are three mandatory fields, one optional-needed field (as in most cases) ...
% and one optional field:
options.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','Continuation of Duffing Equation');
% Mandatory fields: - 'order':  Describes the order of the ODE. As CoSTAR takes a system of differential equations of first order ...
%                               ( dz/dt = f(t,z,param) ), the order is 1.
%                   - 'dim':    Dimension of the state space. In this case, the dimension of the state variable vector z is 2.
%                               -> Allowed values: positive integer (scalar)
%                   - 'rhs':    Right-hand side of equation (2), which we already assigned to the variable Fcn.
%                               -> Allowed values: function handle @(t,z,param)
% Opt-Need. field:  - 'param':  This field expects the previously defined "param" array. We NEED to set it because it is required ...
%                               by Fcn. Moreover, when executing a continuation, param is always necessary. However, param CAN be ...
%                               omitted when only calculating the initial solution (i.e. no continuation is carried out).
%                               -> Allowed values: [1xp] cell array, where p is the number of parameters for the RHS apart from t and z.
% Optional field:   - 'info':   Contains your individual text about the computation. This does not need to be set. However, it may ...
%                               help you later if you have forgotten what this computation is about (default: no text).
%                               -> Allowed values: char or string
%
% Next, there is the "options.opt_sol" structure. There are three mandatory fields and a total of four optional-needed fields. ...
% However, we only need three of the optional-needed fields when computing periodic solutions. The remaining optional-needed ...
% field ('auto_freq') is used instead of 'non_auto_freq' when the system is autonomous. Apart from that, there are two optional ...
% fields, but we use their default values and therefore do not set them.
options.opt_sol = costaropts('sol_type','periodic','approx_method','shooting','cont','on','stability','on', ...
                             'non_auto_freq',non_auto_freq,'act_param',active_parameter);
% Mandatory fields: - 'sol_type':       We need to specify what type of solution CoSTAR has to calculate. In this case, we want to ...
%                                       calculate a periodic solution. It is also possible to use 'ps' instead of 'periodic'.                   
%                   - 'cont':           As we want to execute a continuation, we set 'cont' = 'on'. If 'cont' = 'off', CoSTAR would ...
%                                       only calculate the initial solution z(tau,mu0).
%                   - 'stability':      CoSTAR computes the stability of the solution if 'stability' = 'on'. Set it to 'off' if ...
%                                       it can be skipped. Stability computation is addressed in a separate tutorial in detail.
% Opt-Need. field:  - 'approx_method':  Defines the approximation method used to calculate approximate solutions of equation (2). ...
%                                       For the SHM, the key is 'shooting' or 'shm'.
%                   - 'non_auto_freq':  Function handle setting the non-autonomous frequency of the system. 'non_auto_freq' must ...
%                                       be set if the system exhibits an excitation frequency and is not allowed otherwise.
%                                       -> Allowed values: double or function handle @(mu) with value >= 'freq_limit' (see below)
%                                       -> Default value:  (no default value)
%                   - 'act_param':      As we want to do a continuation, we have to tell CoSTAR where the continuation parameter is ...
%                                       located within param. To do that, we use the previously defined variable "active_parameter". ... 
%                                       Here, this field NEEDS to be set, because we provided a "param" array. If there is no ...
%                                       "param" array, 'act_param' may not be set.
%                                       -> Allowed values: positive integer (scalar) <= length(options.system.param)
%                                       -> Default value:  (no default value)
% Optional fields:  - 'freq_limit':     Defines the bottom limit of allowed base frequency values (frequency limit) of a solution. ...
%                                       CoSTAR stops if the base frequency falls below the limit. 
%                                       -> Allowed values: positive double (scalar)
%                                       -> Default value:  1e-4
%                   - 'display':        Controls the command window output.
%                                       -> Allowed values: 
%                                           - 'off':   Only warnings and errors are displayed.
%                                           - 'final': The latest solution and the termination message are displayed.
%                                           - 'iter':  Displays the iteration process of the initial solution, the latest five ...
%                                                      solutions during continuation and the termination message.
%                                           - 'iter-detailed': Similar to 'iter', but all computed solutions are displayed.
%                                           - 'step-control':  Similar to 'iter-detailed', but with additional step control ...
%                                                              information (if enabled).
%                                           - 'error-control': Similar to 'iter-detailed', but with additional error control ...
%                                                              information (if enabled).
%                                           - 'full':  All possible information are displayed. This is similar to the log file.
%                                       -> Default value: 'iter'
%
% Going on, we have to set the "options.opt_init" structure. The fields of the "options.opt_init" structure depend on the solution ...
% type as well as the chosen approximation method. For periodic solutions using the SHM, there is only one mandatory field and no ...
% optional fields.
options.opt_init = costaropts('ic',IC);
% Mandatory fields: - 'ic':  Defines a point in state space where fsolve starts to iterate towards a point on the periodic orbit ... 
%                            of sought solution z(tau,mu). 'ic' should be as close as possible to the periodic orbit.
%                            -> Allowed values: [dim x 1] (double) array
%                            -> Default value:  (no default value)
%
% So far we have defined all "options" structures which CoSTAR always needs. Since we want to calculate a periodic solution and ...
% execute a continuation, we also have to set the "options.opt_approx_method" as well as the "options.opt_cont" structures.
% The fields of the "options.opt_approx_method" structure depend on the solution type as well as the chosen approximation method. ...
% For periodic solutions using the SHM, there are no mandatory fields and only two optional fields.
options.opt_approx_method = costaropts('solver','ode45','n_shoot',2);
% Optional fields: - 'solver':   Sets the Matlab-proprietary numerical integration algorithm to obtain the time-dependent solution ...
%                                z(tau,mu). Use specialised solvers for stiff ODEs, e.g. ode15s. ...
%                                -> Allowed values: 'ode45', 'ode78', 'ode89', 'ode23', 'ode113', 'ode15s', 'ode23s', 'ode23t', 'ode23tb'
%                                -> Default value:  'ode45'
%                   - 'n_shoot': Defines the number of intervals into which the period of a solution is divided for multiple shooting. ...
%                                This means that n_shoot points on the periodic solution in state space are computed. ...
%                                n_shoot = 1 is denoted as single shooting method, while n_shoot > 1 represents multiple shooting.
%                                -> Allowed values: Positive integer
%                                -> Default value:  2
% NOTE: It is not necessary to set the fields 'solver' and 'n_shoot' since the default values are used. However, it is demonstrated here ...
%       in order to show the fields and possible values.
%
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.05);
% Mandatory fields: - 'mu_limit':  Sets the limits of the continuation. For this purpose, we defined the "mu_limit" variable.
% Optional fields:  - 'step_width':          Defines the initial step width (default: 0.1). The step width will be altered in the ...
%                                            range of step_width*[0.2, 5].
% NOTE: The field belonging to the step control ('step_width') is set to compute the solution curve with a decent resolution. ...
%       Details are explained in a separate step control tutorial so you can ignore this at the moment.
%
% Finally, we are done defining the required settings. All solution type and approximation method specific fields, which are the ...
% fields of "options.opt_init" and "options.opt_approx_method", were explained above. Concerning the rest of the "options" structures ...
% ("options.system", "options.opt_sol" and "options.opt_cont"), only the necessary fields were defined and explained to some extend.
% If you want to have a deeper insight into the "options" structures and its fields, please use the "costarhelp" function by typing ...
% "costarhelp.options" in the command window. In order to directly open the help pages of particular "options" structures, type ...
% "costarhelp.<name_of_options_structure>", e.g. "costarhelp.opt_cont". In case of "options.opt_init" and "options.opt_approx_method", ...
% type "costarhelp.<name_of_options_structure>('PS','SHM')" ("PS": Periodic Solution).


%     1.2.3 Calling CoSTAR and running the simulation
%
% Now we can run the computation. This is done by invoking the "costar" function and passing the "options" structure as an input ...
% argument. "costar" returns two objects: 
% - DynamicalSystem object "DYN":  Saves all the information and settings which are contained in "options". It can be used to ...
%                                  restart the computation and it is necessary for postprocessing.
% - Solution object "S": Stores the information that CoSTAR calculated and can further be used for postprocessing.
%
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options).
%
% During the computation, CoSTAR displays information in the command window:
% - At the beginning, the iteration process of fsolve trying to find the first point on the curve is shown.
% - As soon as fsolve succeeded, CoSTAR announces "Initial solution found!". ...
%   If the stability of a solution is additionally computed, the message is supplemented by the information (stable) or (unstable).
% - After that, the continuation starts and CoSTAR displays "Iter: <XXX> -- mu = <XXX> -- stepwidth = <XXX>" for the latest five ...
%   computed solutions.
%       * "Iter: <XXX>" depicts the number of the computed point on the curve. Please keep in mind that the initial solution is ...
%          counted as Iter = 1.
%       * "mu = <XXX>" shows the value of the continuation parameter mu for the related point on the curve.
%       * "stepwidth = <XXX>" displays the step width that was used to compute the point on the curve.
%       * "(un)stable" reports whether the solution is stable or unstable if the stability was computed.
% - Finally, CoSTAR reports the reason of termination of the continuation.
% 
% Please note that the command window output described above is the default output and it may change if the options.opt_sol field ...
% 'display' is manually set.
% 
% Furthermore, a log file named "CoSTAR_Log_ID_yyyy-mm-dd_hh.mm.ss.txt" (yyyy-mm-dd and hh.mm.ss show the date and time of day ...
% when the computation started) is created in the current directory of MATLAB, giving a more in-depth view of the computation.


%                1.3 Postprocessing                 %
%
% As you have probably noticed, CoSTAR automatically creates a diagram as well. By default, the maximum of the Euclidean norm of the ...
% state space vectors are plotted against the continuation parameter mu. This means that for every computed approximate solution ...
% z(tau,mu), the Euclidean norm of 200 state space vectors is calculated. These 200 state space vectors are distributed equidistantly ...
% in the interval [0, T] (T: periodic time) for one specific mu. After that, the maximum thereof is taken and plotted against mu. ...
% This can be understood as a mapping of the computed curve from R^(n+1)-dimensional space into 2D space.
%
% In the default plot created by CoSTAR, unstable solutions are depicted in red, while stable solutions are shown in blue. ...
% ATTENTION: Please keep in mind that all solutions are plotted in blue if the stability was not computed! ...
% Moreover, CoSTAR also computed bifurcation points, which are shown in the plot ("FB": fold, pitchfork or transcritical bifurcation). ...
% Details are saved to the field 'bifurcation' of the solution object S.
% It is pointed out again that details of the computation of stability and bifurcation points are addressed in a separate tutorial ...
% and are therefore not explained here.
%
% In order to plot z1(tau,mu) = x(tau,mu), we need to call the "solplot" function. Similar to the "costar" function, solplot ...
% expects a structure that defines all required options for the plot.
solplot_options_1 = costaropts('zaxis',@(z) z(:,1),'space','time','mu',[0.5,1,2]);
% Mandatory fields: - 'zaxis':  Defines what is plotted on the vertical axis. "@(z) z(:,1)" plots z1 = x.
%                   - 'space':  Specifies the "domain" of the plot. 'time' plots the time behaviour of the solution.
% Optional fields:  - 'mu':     Defines at which mu-values x(tau,mu) is plotted. If the solution at a specified mu-value is not ...
%                               available, CoSTAR takes the solution that is closest to the desired mu-value (default: all ...
%                               calculated approximate solutions x(tau,mu) are plotted)
%
% Now we can call the "solplot" function. It is a function of the solution class object S, which is why we have to call it by ...
% "S.solplot". Apart from the "opt_solplot" structure, solplot also requires the DynamicalSystem class object DYN. 
% solplot creates the desired plot and additionally returns a structure storing the corresponding z-values, time values and mu-values.
% However, the fields of the output structure depends on the options, so it may change in other cases.
solplot_output_1 = S.solplot(DYN,solplot_options_1);
%
% We do not go into further postprocessing details here as this would exceed the scope of this tutorial.
% Comprehensive explanations of the postprocessing functions can be found in the corresponding tutorials ...
% (Tutorial_Postprocessing_<...>: contplot, solget and solplot) and the associated examples provide short exemplary code. 
% An overview of the available options is given via the costarhelp function by typing "costarhelp.costar" in the command window.


%      1.4 Change of the continuation parameter      %
%
% In the example explained above, the continuation parameter was the excitation frequency eta. How does the "options" structure need ...
% to be modified if we want to vary a different parameter, let's say the coefficient of the nonlinear stiffness kappa?
% Think about it before you continue with this tutorial! If you are not sure which parameters or in what sense parameters need to ...
% be adjusted, please scroll up and read the section "1.2.1 Setting useful variables" again.
%
% First of all, we define new limits for the continuation and set the starting value kappa0:
mu_limit_kappa = [0, 1];
kappa0 = mu_limit_kappa(1);
% We do not change the value of the Damping factor D (= 0.05) and of the excitation amplitude g (= 1). Of course, the function ...
% handle Fcn also remains.
% In contrast to the continuation above, we need to fix the excitation frequency eta since CoSTAR accepts only one continuation ...
% parameter, which is kappa in this case. We set
eta_kappa = 1.5;
% Furthermore, we need to define a new function handle "non_auto_freq_kappa" and assign the specified excitation frequency ...
% eta_kappa to it:
non_auto_freq_kappa = @(mu) eta_kappa;
% Now we can define a new "param" array and set the active parameter:
param_kappa = {kappa0, D, eta_kappa, g};
active_parameter_kappa = 1;
%
% All necessary adjustments to the (important) parameters have been made. We do not need to change the initial condition 'ic' or ...
% the solver. Now we can define the "options" structures and run CoSTAR.
options_kappa.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param_kappa,'info','Continuation of Duffing Equation - kappa');
options_kappa.opt_sol = costaropts('sol_type','periodic','approx_method','shooting','cont','on','stability','on', ...
                                   'non_auto_freq',non_auto_freq_kappa,'act_param',active_parameter_kappa);
options_kappa.opt_init = costaropts('ic',IC);
options_kappa.opt_approx_method = costaropts();     % We use the default solver ode45 so we do not need to define the solver
options_kappa.opt_cont = costaropts('mu_limit',mu_limit_kappa);
%
[S_kappa,DYN_kappa] = costar(options_kappa); 
%
% That's it! Now you also know how to set up CoSTAR if the excitation frequency is not the continuation parameter.


% Finally, we are done with the "Duffing Oscillator" example!
% Going on, we will learn how to deal with autonomous systems on the basis of the van der Pol oscillator.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Ex. No 2: van der Pol Oscillator         %%
%                 Autonomous System                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear variables; clc; close all;                    % Let's clean our "desk" before we begin with our second example.

%                2.1 Introduction                   %
%
% The second example we have a look on is the van der Pol oscillator which is given by the dimensionless differential equation
%   x'' + epsilon * (x^2 - 1) * x' + x = 0   (3)
% epsilon is the coefficient of the non-linearity, x' = dx/dtau and x'' = dx'/dtau (tau being the dimensionless time). However, ...
% CoSTAR always needs the equation to be of the form dz/dt = f(t,z,param) (param represents all additional parameters), which is ...
% why we have to transform (3) into a system of two differential equations of first order. To do that, we introduce the state ...
% variable vector z = [z1; z2] = [x; x'], leading to:
% dz/dtau = f(z,epsilon)  (4)  <=>  z1' = z2
%                                   z2' = - epsilon * (z1^2 - 1) * z2 - z1
%
% Our aim is to approximate the solution of equation (4) z(tau) for epsilon in the range of [0.1, 2.5] using the SHM. As the ...
% system is autonomous, the (angular) frequency of the periodic solution is not known beforehand and has to be calculated as well.
% This entire problem can be rewritten to finding the zero set of the function g: R^(n+2) -> R^(n+1) with g = 0 and for epsilon ...
% in [0.1, 2.5], which defines a curve in (n+2)-dimensional space (n equals the dimension of the system. Here, n = 2 holds). ...
% Hence, a solution z(tau) of (4) is also denoted as a "point on the curve" in the following explanations. ...
% CoSTAR will approximate discrete points on this curve by executing a path continuation.


%       2.2 Continuation of nonlinear damping       %
%
%          2.2.1 Setting useful variables           
%
% Before we specify the settings needed by CoSTAR, we introduce some important variables, which we will use later on. This is not ...
% necessarily needed but it helps to keep an overview of the most important variables/settings. Furthermore, important aspects of ...
% CoSTAR are explained here, making the following section "2.2.2 The CoSTAR settings" clearer.
mu_limit = [0.1, 2.5];                              % mu_limits = epsilon_limits = [0.1, 2.5] sets the limits of our continuation.
%
epsilon0 = mu_limit(1);                             % epsilon0 defines the epsilon-value at which the continuation begins.
Fcn = @(t,z,param) vdP_auto_ap(t,z,param);          % Fcn is a function handle containing the right-hand side (RHS) of equation (4).
% You may have noticed that Fcn is a function of t, z and param, whereas the right-hand side of (4) is only a function of z and ...
% epsilon. This is because CoSTAR ALWAYS requires the definition of Fcn as "Fcn = @(t,z,param) ..." when calculating periodic ...
% solutions. param is a cell array storing all parameters, apart from t and z, that need to be passed to (4). The time t as well ...
% as the state variable vector z are always treated separately and therefore not included in param. Here, we need to define param ...
% since we need to pass the parameter epsilon to Fcn. The name of the "time argument" being t (and not tau) is required by CoSTAR. ...
% It is equal to our dimensionless time tau.
%
% TASK: Open the function "vdP_auto_ap" (click the right mouse button right on it and choose 'Open "vdP_auto_ap"' or go into the ...
%       "RHS" subfolder located in the main path of CoSTAR and double click on the function). Have a look at how the equation is  ...
%       defined and compare that to our system right next to equation (4).
% IMPORTANT: Please pay attention to how the state variables z1 and z2 are defined and how epsilon is extracted from param.
%
% Now we define the "param" array. It is important that it is a cell array AND the order of the parameters placed in param must ...
% correspond to the definitions in the "vdP_auto_ap" function! Since epsilon is the only parameter of the RHS of (4) (apart from z), ...
% we do not need to pay attention to the order of the parameters in param.
param = {epsilon0};                                 % As the continuation starts at epsilon0, we need to place epsilon0 into param
% Next, we have to tell CoSTAR where the continuation parameter is located within param. We also call this the "active parameter". ...
% This is done by defining an integer variable which is used to get the continuation parameter from param by cell indexing, i.e. ...
% continuation parameter mu = param{active_parameter}. As epsilon (epsilon0) is the FIRST (and only) element of param, we set:
active_parameter = 1;
%
% Due to the system being autonomous, we have to set an initial "guess" for the autonomous frequency, which CoSTAR uses as an ...
% initial condition for the autonomous frequency when iterating to the first point on the curve. This parameter is called ...
% "auto_freq" and we set it to
auto_freq = 1;
% Last but not least, we set a variable called "initial condition" (IC). Using the SHM, the initial condition defines a point in ...
% state space where the nonlinear system solver fsolve starts to find the first approximate solution of equation (4) at mu0.
% Remember that the SHM iterates towards a point z_po on the periodic orbit (which is a closed trajectory) of sought solution ...
% z(tau,mu) in state space. The time-dependent approximate solution z(tau,mu) is obtained by time integration starting at z_po.
% The initial condition should be as close as possible to the periodic orbit.
IC = [2; 0];


%            2.2.2 The CoSTAR settings              
%
% All information and settings which CoSTAR needs are defined in a structure array called "options". 
% The structure "options" itself consists of structure arrays that define important parameters used by different modules of CoSTAR.
% 
% CoSTAR ALWAYS needs the following "option" structures: "options.system", "options.opt_sol" and "options.opt_init".
% Depending on the system and problem, CoSTAR might need some or all of the following "option" structures additionally: 
% "options.opt_approx_method":  Required when calculating (quasi-) periodic solutions.
% "options.opt_cont":           Required when running a continuation.
% "options.opt_stability":      Optional when calculating the stability of solutions.
% In this case, "options.opt_approx_method" and "options.opt_cont" are additionally needed.
% 
% Each of the "option" structures contain mandatory fields which ALWAYS have to be set. 
% Furthermore, there are optional fields that CAN be set. Most of the optional fields have a default value. 
% However, even optional fields NEED to be set in some cases (denoted as "Opt-Need.", e.g. 'param' in this example, see below).
% 
% IMPORTANT: All option structures must be created by the "costaropts" function. The syntax of the arguments of "costaropts" is ...
%            equivalent to the MATLAB struct function, i.e. costaropts(fieldname1,value1,...,fieldnameN,valueN).
%
% Let us start with the "options.system" structure. There are three mandatory fields, one optional-needed field (as in most cases) ...
% and one optional field: 
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','Continuation of van der Pol Oscillator');
% Mandatory fields: - 'order':  Describes the order of the ODE. As CoSTAR takes a system of differential equations of first order ...
%                               ( dz/dt = f(t,z,param) ), the order is 1.
%                   - 'dim':    Dimension of the state space. In this case, the dimension of the state variable vector z is 2.
%                               -> Allowed values: positive integer (scalar)
%                   - 'rhs':    Right-hand side of equation (4), which we already assigned to the variable Fcn.   
%                               -> Allowed values: function handle @(t,z,param)
% Opt-Need. field:  - 'param':  This field expects the previously defined "param" array. We NEED to set it because it is required ...
%                               by Fcn. Moreover, when executing a continuation, param is always necessary. However, param CAN be ...
%                               omitted when only calculating the initial solution (i.e. no continuation is carried out).
%                               -> Allowed values: [1xp] cell array, where p is the number of parameters for the RHS apart from t and z.
% Optional field:   - 'info':   Contains your individual text about the computation. This does not need to be set. However, it may ...
%                               help you later if you have forgotten what this computation is about (default: no text).
%                               -> Allowed values: char or string
%
% Next, there is the "options.opt_sol" structure. There are three mandatory fields and a total of four optional-needed fields. ...
% However, we only need three of the optional-needed fields when computing periodic solutions. The remaining optional-needed ...
% field ('non_auto_freq') is used instead of 'auto_freq' when the system is non-autonomous. Apart from that, there are two ...
% optional fields, but we use their default values and therefore do not set them.
options.opt_sol = costaropts('sol_type','periodic','approx_method','shooting','cont','on','stability','on', ...
                             'auto_freq',auto_freq,'act_param',active_parameter);
% Mandatory fields: - 'sol_type':       We need to specify what type of solution CoSTAR has to calculate. In this case, we want to ...
%                                       calculate a periodic solution. It is also possible to use 'ps' instead of 'periodic'.                   
%                   - 'cont':           As we want to execute a continuation, we set 'cont' = 'on'. If 'cont' = 'off', CoSTAR would ...
%                                       only calculate the initial solution z(tau,mu0).
%                   - 'stability':      CoSTAR computes the stability of the solution if 'stability' = 'on'. Set it to 'off' if ...
%                                       it can be skipped. Stability computation is addressed in a separate tutorial in detail.
% Opt-Need. field:  - 'approx_method':  Defines the approximation method used to calculate approximate solutions of equation (4). ...
%                                       For the SHM, the key is 'shooting' or 'shm'.
%                   - 'auto_freq':      Sets the initial value for the autonomous frequency of the system. 'auto_freq' must be ...
%                                       set if the system exhibits an autonomous frequency and is not allowed otherwise.
%                                       -> Allowed values: double (scalar) >= 'freq_limit' (see below)
%                                       -> Default value:  (no default value)
%                   - 'act_param':      As we want to do a continuation, we have to tell CoSTAR where the continuation parameter is ...
%                                       located within param. To do that, we use the previously defined variable "active_parameter". ... 
%                                       Here, this field NEEDS to be set, because we provided a "param" array. If there is no ...
%                                       "param" array, 'act_param' may not be set.
%                                       -> Allowed values: positive integer (scalar) <= length(options.system.param)
%                                       -> Default value:  (no default value)
% Optional fields:  - 'freq_limit':     Defines the bottom limit of allowed base frequency values (frequency limit) of a solution. ...
%                                       CoSTAR stops if the base frequency falls below the limit. 
%                                       -> Allowed values: positive double (scalar)
%                                       -> Default value:  1e-4
%                   - 'display':        Controls the command window output.
%                                       -> Allowed values: 
%                                           - 'off':   Only warnings and errors are displayed.
%                                           - 'final': The latest solution and the termination message are displayed.
%                                           - 'iter':  Displays the iteration process of the initial solution, the latest five ...
%                                                      solutions during continuation and the termination message.
%                                           - 'iter-detailed': Similar to 'iter', but all computed solutions are displayed.
%                                           - 'step-control':  Similar to 'iter-detailed', but with additional step control ...
%                                                              information (if enabled).
%                                           - 'error-control': Similar to 'iter-detailed', but with additional error control ...
%                                                              information (if enabled).
%                                           - 'full':  All possible information are displayed. This is similar to the log file.
%                                       -> Default value: 'iter'
%
% Going on, we have to set the "options.opt_init" structure. The fields of the "options.opt_init" structure depend on the solution ...
% type as well as the chosen approximation method. For periodic solutions using the SHM, there is only one mandatory field and no ...
% optional fields.
options.opt_init = costaropts('ic',IC);
% Mandatory fields: - 'ic':  Defines a point in state space where fsolve starts to iterate towards a point on the periodic orbit ... 
%                            of sought solution z(tau,mu). 'ic' should be as close as possible to the periodic orbit.
%                            -> Allowed values: [dim x 1] (double) array
%                            -> Default value:  (no default value)
%
% So far we have defined all "options" structures which CoSTAR always needs. Since we want to calculate a periodic solution and ...
% execute a continuation, we also have to set the "options.opt_approx_method" as well as the "options.opt_cont" structures.
% The fields of the "options.opt_approx_method" structure depend on the solution type as well as the chosen approximation method. ...
% For periodic solutions using the SHM, there are no mandatory fields and only two optional fields.
options.opt_approx_method = costaropts('solver','ode45','n_shoot',2);
% Optional fields: - 'solver':   Sets the Matlab-proprietary numerical integration algorithm to obtain the time-dependent solution ...
%                                z(tau,mu). Use specialised solvers for stiff ODEs, e.g. ode15s. ...
%                                -> Allowed values: 'ode45', 'ode78', 'ode89', 'ode23', 'ode113', 'ode15s', 'ode23s', 'ode23t', 'ode23tb'
%                                -> Default value:  'ode45'
%                   - 'n_shoot': Defines the number of intervals into which the period of a solution is divided for multiple shooting. ...
%                                This means that n_shoot points on the periodic solution in state space are computed. ...
%                                n_shoot = 1 is denoted as single shooting method, while n_shoot > 1 represents multiple shooting.
%                                -> Allowed values: Positive integer
%                                -> Default value:  2
% NOTE: It is not necessary to set the fields 'solver' and 'n_shoot' since the default values are used. However, it is demonstrated here ...
%       in order to show the fields and possible values.
%
options.opt_cont = costaropts('mu_limit',mu_limit);
% Mandatory fields: - 'mu_limit':  Sets the limits of the continuation. For this purpose, we defined the "mu_limit" variable.
%
% Finally, we are done defining the required settings. All solution type and approximation method specific fields, which are the ...
% fields of "options.opt_init" and "options.opt_approx_method", were explained above. Concerning the rest of the "options" structures ...
% ("options.system", "options.opt_sol" and "options.opt_cont"), only the necessary fields were defined and explained to some extend.
% If you want to have a deeper insight into the "options" structures and its fields, please use the "costarhelp" function by typing ...
% "costarhelp.options" in the command window. In order to directly open the help pages of particular "options" structures, type ...
% "costarhelp.<name_of_options_structure>", e.g. "costarhelp.opt_cont". In case of "options.opt_init" and "options.opt_approx_method", ...
% type "costarhelp.<name_of_options_structure>('PS','SHM')" ("PS": Periodic Solution).


%     2.2.3 Calling CoSTAR and running the simulation
%
% Now we can run the computation. This is done by invoking the "costar" function and passing the "options" structure as an input ...
% argument. "costar" returns two objects: 
% - DynamicalSystem object "DYN":  Saves all the information and settings which are contained in "options". It can be used to ...
%                                  restart the computation and it is necessary for postprocessing.
% - Solution object "S": Stores the information that CoSTAR calculated and can further be used for postprocessing.
%
[S,DYN] = costar(options);                          % Calling CoSTAR and performing the continuation
%
% During the computation, CoSTAR displays information in the command window:
% - At the beginning, the iteration process of fsolve trying to find the first point on the curve is shown.
% - As soon as fsolve succeeded, CoSTAR announces "Initial solution found!". ...
%   If the stability of a solution is additionally computed, the message is supplemented by the information (stable) or (unstable).
% - After that, the continuation starts and CoSTAR displays "Iter: <XXX> -- mu = <XXX> -- stepwidth = <XXX>" for the latest five ...
%   computed solutions.
%       * "Iter: <XXX>" depicts the number of the computed point on the curve. Please keep in mind that the initial solution is ...
%          counted as Iter = 1.
%       * "mu = <XXX>" shows the value of the continuation parameter mu for the related point on the curve.
%       * "stepwidth = <XXX>" displays the step width that was used to compute the point on the curve.
%       * "(un)stable" reports whether the solution is stable or unstable if the stability was computed.
% - Finally, CoSTAR reports the reason of termination of the continuation.
% 
% Please note that the command window output described above is the default output and it may change if the options.opt_sol field ...
% 'display' is manually set.
% 
% Furthermore, a log file named "CoSTAR_Log_ID_yyyy-mm-dd_hh.mm.ss.txt" (yyyy-mm-dd and hh.mm.ss show the date and time of day ...
% when the computation started) is created in the current directory of MATLAB, giving a more in-depth view of the computation.


%                2.3 Postprocessing                 %
%
% As you have probably noticed, CoSTAR automatically creates a diagram as well. By default, the maximum of the Euclidean norm of the ...
% state space vectors are plotted against the continuation parameter mu. This means that for every computed approximate solution ...
% z(tau,mu), the Euclidean norm of 200 state space vectors is calculated. These 200 state space vectors are distributed equidistantly ...
% in the interval [0, T] (T: periodic time) for one specific mu. After that, the maximum thereof is taken and plotted against mu. ...
% This can be understood as a mapping of the computed curve from R^(n+2)-dimensional space into 2D space.
%
% In the default plot created by CoSTAR, unstable solutions are depicted in red, while stable solutions are shown in blue. ...
% ATTENTION: Please keep in mind that all solutions are plotted in blue if the stability was not computed! ...
% Moreover, CoSTAR can also compute bifurcation points, but there are no bifurcation points occurring in this example. ...
% It is pointed out again that details of the computation of stability and bifurcation points are addressed in a separate tutorial ...
% and are therefore not explained here.
%
% In order to plot z1(tau,mu) = x(tau,mu), we need to call the "solplot" method. Similar to the "costar" function, solplot expects a ...
% structure that defines all required options for the plot.
solplot_options_2 = costaropts('zaxis',@(z) z(:,1),'space','time','mu',[1,2]);
% Mandatory fields: - 'zaxis':  Defines what is plotted on the vertical axis. "@(z) z(:,1)" plots z1 = x.
%                   - 'space':  Specifies the "domain" of the plot. 'time' plots the time behaviour of the solution.
% Optional fields:  - 'mu':     Defines at which mu-values x(tau,mu) is plotted. If the solution at a specified mu-value is not ...
%                               available, CoSTAR takes the solution that is closest to the desired mu-value (default: all ...
%                               calculated approximate solutions x(tau,mu) are plotted).
%
% Now we can call the "solplot" method. It is a method of the solution class object S, which is why we have to call it by ...
% "S.solplot". Apart from the "opt_solplot" structure, solplot also requires the DynamicalSystem class object DYN. 
% solplot creates the desired plot and additionally returns a structure storing the corresponding z-values, time values and mu-values.
% However, the fields of the output structure depends on the options, so it may change in other cases.
solplot_output_2 = S.solplot(DYN,solplot_options_2);
%
% We do not go into further postprocessing details here as this would exceed the scope of this tutorial.
% Comprehensive explanations of the postprocessing functions can be found in the corresponding tutorials ...
% (Tutorial_Postprocessing_<...>: contplot, solget and solplot) and the associated examples provide short exemplary code. 
% An overview of the available options is given via the costarhelp function by typing "costarhelp.costar" in the command window.


% Finally, we are done with the "van der Pol Oscillator" example!



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Final Words                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The CoSTAR tutorial on calculating periodic solutions using the (Multiple) Shooting Method is now finished.
% For additional information, please use the "costarhelp" function and/or the CoSTAR manual.

% If you are interested in learning about further capabilities of CoSTAR, you are invited to have a look at the other tutorials as well.

% See you soon!