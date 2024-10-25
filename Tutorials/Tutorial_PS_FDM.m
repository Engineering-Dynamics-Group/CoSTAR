%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      CoSTAR                       %
%   Continuation of Solution Torus AppRoximations   %
%                                                   %
%                    Tutorial:                      %
%               Periodic Solutions                  %
%          - Finite Difference Method -             %
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
% In this tutorial, we will look at how to compute periodic solutions using the Finite Difference Method (FDM) to approximate the ...           
% solution (associated example: Example_PS_FDM). We will learn how to set up CoSTAR and how to correctly define all the required ...
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
%       ->  FDM specific: - Usage of the "options.opt_init"-fields:
%                               * 'c1' and 's1' in section 1.2
%                               * 'fdm_sol' in section 1.4
%                         - Usage of the "options.opt_approx_method"-fields:
%                               * 'n_int', 'scheme' and 'approx_order' in section 1.2
%                               * 'n_int' and 'points' in section 1.4
%
%  2. Continuation of the van der Pol oscillator (autonomous system)
%       2.1 Introduction
%       2.2 Continuation of nonlinear damping
%           2.2.1 Setting useful variables
%           2.2.2 The CoSTAR settings
%           2.2.3 Calling CoSTAR and running the simulation
%       2.3 Postprocessing: Plot of time-dependent behaviour of solution x(tau,mu)
%       ->  FDM specific: - Usage of the "options.opt_init"-fields 'c1' and 's1'
%                         - Usage of the "options.opt_approx_method"-fields 'n_int' and 'points'
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
% Our aim is to approximate the solution of equation (2) z(tau) for eta in the range of [0.01, 2.5] using the FDM. This problem ...
% can be rewritten to finding the zero set of the function g: R^(n+1) -> R^n with g = 0 and for eta in [0.01, 2.5], which ...
% defines a curve in (n+1)-dimensional space (n depends on the dimension of the system as well as the finite difference discretisation). 
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
% parameter. In order to deal with this uncertainty, we introduce a function handle which sets the non-autonomous frequency. In this ...
% case, the non-autonomous frequency is the continuation parameter mu, so we set
non_auto_freq = @(mu) mu;
% "non_auto_freq" ALWAYS have to be a function handle @(mu), even if it is not dependent on mu! For example, if we want to vary the ...
% damping factor D, the non-autonomous frequency would have to be fixed, let us say eta = 1.5. In that case, we would have to define ...
% the function handle as non_auto_freq = @(mu) 1.5. This rather unconventional approach not only keeps the code as compact as ...
% possible but also allows a great flexibility for the user since the non-autonomous frequency can be an arbitrary function of mu ...
% (e.g. we can set non_auto_freq = @(mu) 10*mu etc.).
%
% Last but not least, we set variables which are used to create an initial value (IV) where the nonlinear system solver fsolve ...
% starts to find the first approximate solution of equation (2) at mu0 (i.e. the first point on the curve), referred to as ...
% z0 = z0(tau,mu0). Using finite differences, the initial value is created via a first order Fourier series:
% z_IV (tau,eta0) = C0 + C1*cos(eta0*tau) + S1*sin(eta0*tau). 
% The coefficients C0, C1 and S1 can be set by the user. They have to be column vectors of dimension dim and their default value is ...
% the zero vector. It is important to note that the better z_IV matches the approximate solution z0, the faster fsolve converges ...
% towards z0. If z_IV is too far away from z0 and/or the discretisation is too coarse, it may happen that fsolve does not converge ...
% at all! For our Duffing example, we use z1_IV (tau,eta0) = x_IV (tau,eta0) = g*cos(eta0*tau). Consequently, ...
% z2_IV (tau,eta0) = x'_IV (tau,eta0) = - eta0*g*sin(eta0*tau). Thus, we need to define C1 and S1 as follows:
C1 = [g; 0];
S1 = [0; -g*eta0];


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
% Let us start with the "options.system" structure:
options.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','Continuation of Duffing Equation');
% Mandatory fields: - 'order':  Describes the order of the ODE. As CoSTAR takes a system of differential equations of first order ...
%                               ( dz/dt = f(t,z,param) ), the order is 1.
%                   - 'dim':    Dimension of the state space. In this case, the dimension of the state variable vector z is 2.
%                   - 'rhs':    Right-hand side of equation (2), which we already assigned to the variable Fcn.                
% Opt-Need. field:  - 'param':  This field expects the previously defined "param" array. We NEED to set it because it is required ...
%                               by Fcn. Moreover, when executing a continuation, param is always necessary. However, param CAN be ...
%                               omitted when only calculating the initial solution (i.e. no continuation is carried out).
% Optional field:   - 'info':   Contains your individual text about the computation. This does not need to be set. However, it may ...
%                               help you later if you have forgotten what this computation is about (default: no text).
%
% Next, there is the "options.opt_sol" structure:
options.opt_sol = costaropts('sol_type','periodic','approx_method','finite-difference','cont','on','stability','on', ...
                             'non_auto_freq',non_auto_freq,'act_param',active_parameter);
% Mandatory fields: - 'sol_type':       We need to specify what type of solution CoSTAR has to calculate. In this case, we want to ...
%                                       calculate a periodic solution. It is also possible to use 'ps' instead of 'periodic'.
%                   - 'approx_method':  Defines the approximation method used to calculate approximate solutions of equation (2). ...
%                                       For the FDM, the key is 'finite-difference' or 'fdm'.
%                   - 'cont':           As we want to execute a continuation, we set 'cont' = 'on'. If 'cont' = 'off', CoSTAR would ...
%                                       only calculate the initial solution z(tau,mu0).
%                   - 'stability':      CoSTAR computes the stability of the solution if 'stability' = 'on'. Stability is addressed ...
%                                       in a separate tutorial in detail.
% Opt-Need. field:  - 'non_auto_freq':  Function handle setting the non-autonomous frequency of the system. 'non_auto_freq' must ...
%                                       be set if the system exhibit an excitation frequency and is not allowed otherwise.
%                   - 'act_param':      As we want to do a continuation, we have to tell CoSTAR where the continuation parameter is ...
%                                       located within param. To do that, we use the previously defined variable "active_parameter". ... 
%                                       Here, this field NEEDS to be set, because we provided a "param" array. If there is no ...
%                                       "param" array, 'act_param' may not be set.
%
% Going on, we have to set the "options.opt_init" structure. The fields of the "options.opt_init" structure depend on the solution ...
% type as well as the chosen approximation method. For periodic solutions using the FDM, there are no mandatory fields and four ...
% optional fields of which a maximum of 3 may be used. It is not allowed to use 'fdm_sol' in combination with 'c0', 'c1' or 's1'.
options.opt_init = costaropts('c1',C1,'s1',S1);
% Optional fields: - 'c0':        Sets the Fourier series coefficient of the constant term used to create an initial value for fsolve ...
%                                 to find the first approximate solution z0(tau,mu0). We do not need 'c0' here so we do not specify it.
%                                 -> Allowed values: [dim x 1] (double) array
%                                 -> Default value:  zeros(dim,1)
%                  - 'c1':        Sets the Fourier series coefficient of the "cos(eta0*tau)" term used to create an initial value for ...
%                                 fsolve to find the first approximate solution z0(tau,mu0). We already defined C1 above.
%                                 -> Allowed values: [dim x 1] (double) array
%                                 -> Default value:  zeros(dim,1)
%                  - 's1':        Sets the Fourier series coefficient of the "sin(eta0*tau)" term used to create an initial value for ...
%                                 fsolve to find the first approximate solution z0(tau,mu0). We already defined S1 above.
%                                 -> Allowed values: [dim x 1] (double) array
%                                 -> Default value:  zeros(dim,1)
%                   - 'fdm_sol':  The initial value can be taken from an already calculated solution. 'fdm_sol' takes the method ...
%                                 solution vector s (stored in Solution_object.s). If the number of intervals 'n_int', defined in ...
%                                 opt_approx_method (or its default value, see below), does not match the number of intervals of ...
%                                 'fdm_sol' (n_int_fdm_sol), the provided solution is interpolated. Moreover, the fields 'c0', 'c1' ...
%                                 and 's1' are not allowed when providing 'fdm_sol'.
%                                 -> Allowed values [dim*n_int_fdm_sol x 1] (double) array
%                                 -> Default value: (no default value)
% NOTE: For usage of the field 'fdm_sol', please see section 1.4.
%                  
% So far we have defined all "options" structures which CoSTAR always needs. Since we want to calculate a periodic solution and ...
% execute a continuation, we also have to set the "options.opt_approx_method" as well as the "options.opt_cont" structures.
% The fields of the "options.opt_approx_method" structure depend on the solution type as well as the chosen approximation method. ...
% For periodic solutions using the FDM, there are no mandatory fields and four optional field of which a maximum of three may be used.
options.opt_approx_method = costaropts('n_int',100,'scheme','central','approx_order',6);
% Optional fields: - 'n_int':         Sets the number of intervals into which one period is divided. That means that the solution ...
%                                     z(tau,mu) is approximated by (n_int+1) points in the range of [0, T] (T: periodic time). ... 
%                                     This defines the resolution of the finite difference grid.
%                                     -> Allowed values: Positive integer >= 2
%                                     -> Default value:  100
%                  - 'scheme':        This field refers to the discretisation scheme used to approximate the derivation dz/dtau. ...
%                                     -> Allowed values: 'central', 'forward', 'backward'
%                                     -> Default value:  'central'
%                  - 'approx_order':  Defines the order of approximation of dz/dtau.
%                                     -> Allowed values: Positive integer. Must be even number if 'scheme' = 'central'
%                                     -> Default value:  6
%                  - 'points':        Specifies the local grid points which are used to approximate dz/dtau. E.g.: If 'points' is ...
%                                     set to [-1,0,1], dz/dtau at grid point index i is approximated by the state space vectors at ...
%                                     (i-1), (i+0) and (i+1). Can be used INSTEAD OF 'scheme' and 'approx_order', i.e. 'points' is ...
%                                     not allowed when specifying 'scheme' and/or 'approx_order' and vice versa.
%                                     -> Allowed values: Integer vector
%                                     -> Default value:  [-3,-2,-1,0,1,2,3]                                                                        
% NOTE: It is not necessary to set 'n_int', 'scheme' and 'approx_order' since the default values are used. However, it is ...
%       demonstrated here in order to show the fields and possible values. For usage of the field 'points', please see section ...
%       1.4 or section 2.2 of the "van der Pol Oscillator" example below.
%
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.25,'step_control_param',[3,7.5/180*pi]);
% Mandatory fields: - 'mu_limit':            Sets the limits of the continuation. For this purpose, we defined the "mu_limit" variable.
% Optional fields:  - 'step_width':          Defines the initial step width (default: 0.1). The step width will be altered in the ...
%                                            range of step_width*[0.2, 5].
%                   - 'step_control_param':  Sets properties for the step control (default value for step control method 'angle', ...
%                                            which is used here: [3,3/180*pi]).
% NOTE: The fields belonging to the step control ('step_width' and 'step_control_param') are set in order to reduce the numerical ...
%       effort. They are explained in a separate step control tutorial so you can ignore them at the moment.
%
% Finally, we are done defining the required settings. All solution type and approximation method specific fields, which are the ...
% fields of "options.opt_init" and "options.opt_approx_method", were explained above. Concerning the rest of the "options" structures ...
% ("options.system", "options.opt_sol" and "options.opt_cont"), only the necessary fields were defined and explained to some extend.
% If you want to have a deeper insight into the "options" structures and its fields, please use the "costarhelp" function by typing ...
% "costarhelp.options" in the command window. In order to directly open the help pages of particular "options" structures, type ...
% "costarhelp.<name_of_options_structure>", e.g. "costarhelp.opt_cont". In case of "options.opt_init" and "options.opt_approx_method", ...
% type "costarhelp.<name_of_options_structure>('PS','FDM')" ("PS": Periodic Solution).


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
% - As soon as fsolve succeeded, CoSTAR announces "Initial solution found!".
% - After that, CoSTAR displays "Iter: <XXX> -- mu = <XXX> -- stepwidth = <XXX>" when a new point on the curve has been calculated.
%       * "Iter" depicts the number of points on the curve which already have been calculated.
%       * "mu" shows the value of the continuation parameter mu at the latest point on the curve.
%       * "stepwidth" displays the step width that was used to calculate the latest point on the curve.
% - When the step width was adapted, CoSTAR shows the new step width and some additional information. You can ignore it at this point.
% - Finally, CoSTAR reports the reason of termination of the continuation.
%
% Apart from the information in the command window, CoSTAR plots the maximum of the Euclidean norm of the state space vectors ...
% against the continuation parameter mu. This means that for every computed approximate solution z(tau,mu), the Euclidean norm ...
% of 200 state space vectors is calculated. These 200 state space vectors are equidistant distributed with respect to time in ...
% the interval [0, T] for one specific mu. After the Euclidean norm of all 200 state space vectors were computed, the maximum ...
% thereof is taken and plotted against mu. This is essentially a projection of the curve in R^(n+1)-dimensional space into 2D space.
%
% Moreover, unstable solutions are depicted in red since we set the field 'stability' of the "options.opt_sol" structure to 'on'. ...
% CoSTAR also computed bifurcation points, which are shown in the plot ("FB": fold, pitchfork or transcritical bifurcation). ...
% Details are saved to the field 'bifurcation' of the solution object S.
% It is pointed out again that details of the computation of stability and bifurcation points are addressed in a separate tutorial ...
% and are therefore not explained here. 


%                1.3 Postprocessing                 %
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
% Instead of the Fourier series coefficients 'c1' and 's1', we use 'fdm_sol' this time. For 'fdm_sol', we take the first ...
% solution point of the continuation above, which is the first column of the field 's' of the solution object S.
% Furthermore, we set 'n_int' = 150, so FDM_sol is interpolated to adapt it to the new number of intervals.
FDM_sol = S.s(:,1);
%
% All necessary adjustments to the (important) parameters have been made. We can now define the "options" structures and run CoSTAR.
options_kappa.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param_kappa,'info','Continuation of Duffing Equation - kappa');
options_kappa.opt_sol = costaropts('sol_type','periodic','approx_method','finite-difference','cont','on','stability','on', ...
                                   'non_auto_freq',non_auto_freq_kappa,'act_param',active_parameter_kappa);
options_kappa.opt_init = costaropts('fdm_sol',FDM_sol);
options_kappa.opt_approx_method = costaropts('n_int',150,'points',[-4,-3,-2,-1,0,1,2]);     % The field 'points' is used this time
options_kappa.opt_cont = costaropts('mu_limit',mu_limit_kappa,'step_width',0.25,'step_control_param',[3,7.5/180*pi]);
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
% Our aim is to approximate the solution of equation (4) z(tau) for epsilon in the range of [0.1, 2.5] using the FDM. As the ...
% system is autonomous, the (angular) frequency of the periodic solution is not known beforehand and has to be calculated as well.
% This entire problem can be rewritten to finding the zero set of the function g: R^(n+2) -> R^(n+1) with g = 0 and for epsilon ...
% in [0.1, 2.5], which defines a curve in (n+2)-dimensional space (n depends on the dimension of the system as well as the finite ...
% difference discretisation). Hence, a solution z(tau) of (4) is also denoted as a "point on the curve" in the following explanations. ...
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
% initial value (IV) for the autonomous frequency when iterating to the first point on the curve. This parameter is called ...
% "auto_freq" and we set it to
auto_freq = 1;
% Last but not least, we set variables which are needed to completely create an initial value where the nonlinear system ...
% solver fsolve starts to find the first approximate solution of equation (4) at mu0 (i.e. the first point on the curve), referred ...
% to as z0 = z0(tau,mu0). Using finite differences, the initial value is created via a first order Fourier series:
% z_IV (tau) = C0 + C1 * cos(auto_freq * tau) + S1 * sin(auto_freq * tau). 
% The coefficients C0, C1 and S1 can be set by the user. They have to be column vectors of dimension dim and their default value is ...
% the zero vector. It is important to note that the better z_IV matches the approximate solution z0, the faster fsolve converges ...
% towards z0. If z_IV is too far away from z0 and/or the discretisation is too coarse, it may happen that fsolve does not converge ...
% at all! For our van der Pol example, we use z1_IV (tau) = x_IV (tau) = 2 * cos(auto_freq * tau). Consequently, ...
% z2_IV (tau) = x'_IV (tau) = - auto_freq * 2 * sin(auto_freq * tau). Thus, we need to define C1 and S1 as follows:
C1 = [2; 0];
S1 = [0; -2];


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
% Let us start with the "options.system" structure: 
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','Continuation of van der Pol Oscillator');
% Mandatory fields: - 'order':  Describes the order of the ODE. As CoSTAR takes a system of differential equations of first order ...
%                               ( dz/dt = f(t,z,param) ), the order is 1.
%                   - 'dim':    Dimension of the state space. In this case, the dimension of the state variable vector z is 2.
%                   - 'rhs':    Right-hand side of equation (4), which we already assigned to the variable Fcn.                
% Opt-Need. field:  - 'param':  This field expects the previously defined "param" array. We NEED to set it because it is required ...
%                               by Fcn. Moreover, when executing a continuation, param is always necessary. However, param CAN be ...
%                               omitted when only calculating the initial solution (i.e. no continuation is carried out).
% Optional field:   - 'info':   Contains your individual text about the computation. This does not need to be set. However, it may ...
%                               help you later if you have forgotten what this computation is about (default: no text).
%
% Next, there is the "options.opt_sol" structure:
options.opt_sol = costaropts('sol_type','periodic','approx_method','finite-difference','cont','on','stability','off', ...
                             'auto_freq',auto_freq,'act_param',active_parameter);
% Mandatory fields: - 'sol_type':       We need to specify what type of solution CoSTAR has to calculate. In this case, we want to ...
%                                       calculate a periodic solution. It is also possible to use 'ps' instead of 'periodic'.
%                   - 'approx_method':  Defines the approximation method used to calculate approximate solutions of equation (4). ...
%                                       For the FDM, the key is 'finite-difference' or 'fdm'.
%                   - 'cont':           As we want to execute a continuation, we set 'cont' = 'on'. If 'cont' = 'off', CoSTAR would ...
%                                       only calculate the initial solution z(tau,mu0).
%                   - 'stability':      CoSTAR computes the stability of the solution if 'stability' = 'on' (this is addressed in a ...
%                                       separate tutorial in detail). We do not need to calculate this, so we set the field to 'off'.
% Opt-Need. field:  - 'auto_freq':      Sets the initial value for the autonomous frequency of the system. 'auto_freq' must be ...
%                                       set if the system exhibits an autonomous frequency and is not allowed otherwise.
%                   - 'act_param':      As we want to do a continuation, we have to tell CoSTAR where the continuation parameter is ...
%                                       located within param. To do that, we use the previously defined variable "active_parameter". ... 
%                                       Here, this field NEEDS to be set, because we provided a "param" array. If there is no ...
%                                       "param" array, 'act_param' may not be set.
%
% Going on, we have to set the "options.opt_init" structure. The fields of the "options.opt_init" structure depend on the solution ...
% type as well as the chosen approximation method. For periodic solutions using the FDM, there are no mandatory fields and four ...
% optional fields of which a maximum of 3 may be used. It is not allowed to use 'fdm_sol' in combination with 'c0', 'c1' or 's1'.
options.opt_init = costaropts('c1',C1,'s1',S1);
% Optional fields: - 'c0':        Sets the Fourier series coefficient of the constant term used to create an initial value for fsolve ...
%                                 to find the first approximate solution z0(tau,mu0). We do not need 'c0' here so we do not specify it.
%                                 -> Allowed values: [dim x 1] (double) array
%                                 -> Default value:  zeros(dim,1)
%                  - 'c1':        Sets the Fourier series coefficient of the "cos(eta0*tau)" term used to create an initial value for ...
%                                 fsolve to find the first approximate solution z0(tau,mu0). We already defined C1 above.
%                                 -> Allowed values: [dim x 1] (double) array
%                                 -> Default value:  zeros(dim,1)
%                  - 's1':        Sets the Fourier series coefficient of the "sin(eta0*tau)" term used to create an initial value for ...
%                                 fsolve to find the first approximate solution z0(tau,mu0). We already defined S1 above.
%                                 -> Allowed values: [dim x 1] (double) array
%                                 -> Default value:  zeros(dim,1)
%                   - 'fdm_sol':  The initial value can be taken from an already calculated solution. 'fdm_sol' takes the method ...
%                                 solution vector s (stored in Solution_object.s). If the number of intervals 'n_int', defined in ...
%                                 opt_approx_method (or its default value, see below), does not match the number of intervals of ...
%                                 'fdm_sol' (n_int_fdm_sol), the provided solution is interpolated. Moreover, the fields 'c0', 'c1' ...
%                                 and 's1' are not allowed when providing 'fdm_sol'.
%                                 -> Allowed values [dim*n_int_fdm_sol x 1] (double) array
%                                 -> Default value: (no default value)
% NOTE: For usage of the field 'fdm_sol', please see section 1.4 of the "Duffing Oscillator" example above.
%  
% So far we have defined all "options" structures which CoSTAR always needs. Since we want to calculate a periodic solution and ...
% execute a continuation, we also have to set the "options.opt_approx_method" as well as the "options.opt_cont" structures.
% The fields of the "options.opt_approx_method" structure depend on the solution type as well as the chosen approximation method. ...
% For periodic solutions using the FDM, there are no mandatory fields and four optional field of which a maximum of three may be used.
options.opt_approx_method = costaropts('n_int',100,'points',[-4,-3,-2,-1,0,1,2]);
% Optional fields: - 'n_int':         Sets the number of intervals into which one period is divided. That means that the solution ...
%                                     z(tau,mu) is approximated by (n_int+1) points in the range of [0, T] (T: periodic time). ... 
%                                     This defines the resolution of the finite difference grid.
%                                     -> Allowed values: Positive integer >= 2
%                                     -> Default value:  100
%                  - 'points':        Specifies the local grid points which are used to approximate dz/dtau. E.g.: If 'points' is ...
%                                     set to [-1,0,1], dz/dtau at grid point index i is approximated by the state space vectors at ...
%                                     (i-1), (i+0) and (i+1).
%                                     -> Allowed values: Integer vector
%                                     -> Default value:  [-3,-2,-1,0,1,2,3]   
%                  - 'scheme':        This field refers to the discretisation scheme used to approximate the derivation dz/dtau. ...
%                                     Can be used INSTEAD OF 'points', i.e. 'scheme' is not allowed when specifying 'points' and ...
%                                     vice versa.
%                                     -> Allowed values: 'central', 'forward', 'backward'
%                                     -> Default value:  'central'
%                  - 'approx_order':  Defines the order of approximation of dz/dtau. Can be used INSTEAD OF 'points', i.e. ...
%                                     'approx_order' is not allowed when specifying 'points' and vice versa.
%                                     -> Allowed values: Positive integer. Must be even number if 'scheme' = 'central'
%                                     -> Default value:  6
% NOTE: It is not necessary to set the field 'n_int' since the default value is used. However, it is demonstrated here in order ...
%       to show the field and a possible value. For usage of the fields 'scheme' and 'approx_order', please see the ...
%       "Duffing Oscillator" example above.
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
% type "costarhelp.<name_of_options_structure>('PS','FDM')" ("PS": Periodic Solution).


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
% - As soon as fsolve succeeded, CoSTAR announces "Initial solution found!".
% - After that, CoSTAR displays "Iter: <XXX> -- mu = <XXX> -- stepwidth = <XXX>" when a new point on the curve has been calculated.
%       * "Iter" depicts the number of points on the curve which already have been calculated.
%       * "mu" shows the value of the continuation parameter mu at the latest point on the curve.
%       * "stepwidth" displays the step width that was used to calculate the latest point on the curve.
% - When the step width was adapted, CoSTAR shows the new step width and some additional information. You can ignore it at this point.
% - Finally, CoSTAR reports the reason of termination of the continuation.
%
% Apart from the information in the command window, CoSTAR plots the maximum of the Euclidean norm of the state space vectors ...
% against the continuation parameter mu. This means that for every computed approximate solution z(tau,mu), the Euclidean norm ...
% of 200 state space vectors is calculated. These 200 state space vectors are equidistant distributed with respect to time in ...
% the interval [0, T] for one specific mu. After the Euclidean norm of all 200 state space vectors were computed, the maximum ...
% thereof is taken and plotted against mu. This is essentially a projection of the curve in R^(n+2)-dimensional space into 2D space.


%               2.3  Postprocessing                 %
%
% In order to plot z1(tau,mu) = x(tau,mu), we need to call the "solplot". Similar to the "costar" function, solplot expects a ...
% structure that defines all required options for the plot.
solplot_options_2 = costaropts('zaxis',@(z) z(:,1),'space','time','mu',[1,2]);
% Mandatory fields: - 'zaxis':  Defines what is plotted on the vertical axis. "@(z) z(:,1)" plots z1 = x.
%                   - 'space':  Specifies the "domain" of the plot. 'time' plots the time behaviour of the solution.
% Optional fields:  - 'mu':     Defines at which mu-values x(tau,mu) is plotted. If the solution at a specified mu-value is not ...
%                               available, CoSTAR takes the solution that is closest to the desired mu-value (default: all ...
%                               calculated approximate solutions x(tau,mu) are plotted).
%
% Now we can call the "solplot" function. It is a function of the solution class object S, which is why we have to call it by ...
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

% The CoSTAR tutorial on calculating periodic solutions using the Finite Difference Method is now finished.
% For additional information, please use the "costarhelp" function and/or the CoSTAR manual.

% If you are interested in learning about further capabilities of CoSTAR, you are invited to have a look at the other tutorials as well.

% See you soon!