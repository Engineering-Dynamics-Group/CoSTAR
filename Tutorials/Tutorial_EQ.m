%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      CoSTAR                       %
%   Continuation of Solution Torus AppRoximations   %
%                                                   %
%                    Tutorial:                      %
%              Equilibrium Solutions                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%           Welcome to the CoSTAR Tutorials!        %
%
% In this tutorial, we will look at how to compute equilibrium solutions using CoSTAR on the basis of two different examples:

%  1. The "Parable" example
%       1.1 Introduction
%       1.2 Continuation 
%           1.2.1 Setting useful variables
%           1.2.2 The CoSTAR settings
%           1.2.3 Calling CoSTAR and running the simulation
%       1.3 Postprocessing: Plot of the parable in a (mu,z)-diagram
%
%  2. The "Pitchfork Bifurcation" example
%       2.1 Introduction
%       2.2 Continuation - Part 1: The solution branch z = 0
%           2.2.1 Setting useful variables
%           2.2.2 The CoSTAR settings
%           2.2.3 Calling CoSTAR and running the simulation
%       2.3 Continuation - Part 2: The top and bottom solution branches for mu >= 0
%       2.4 Postprocessing: Plot of multiple curves / solution branches in a single (mu,z)-diagram

% NOTE: In order to learn how to compute equilibrium solutions, the "Parable" example is sufficient. 
%       The "pitchfork" example extends the "Parable" example by explaining how to deal with multiple solution branches. It ...
%       assumes you have gone through the "Parable" example.
%
% IMPORTANT: Explanations can be found at the beginning of new lines as well as right next to code, so make sure that you do not ...
%            miss some helpful comments!
%            Furthermore, it is advised to execute the code line by line or in blocks when going through this tutorial and do not ...
%            execute the complete script by clicking "Run" or pressing "F5". To do that, you can select the desired lines of code ...
%            to be executed and press "F9" or you can click right on the selected code and choose "Evaluate Selection in Command ...
%            Window". This should help to better understand what is going on and how CoSTAR works. 
%
% Alright, let's clean our "desk" and let's get startet!
clear variables; clc; close all;                                        % clear workspace; clear command window; close all figures

% We need to add the main CoSTAR folder (and all subfolders) to MATLAB's search path. 
% Assuming you are running this script within the "Tutorials" subfolder of CoSTAR, this is done by:
% addpath(genpath('..\'))                                               % genpath() generates a search path containing all subfolders
% (The command is set as comment in order to prevent unwanted folders to be added to MATLAB's search path.)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Example No. 1: Parable             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                1.1 Introduction                   %
%
% The first example we have a look on is called "Parable". As the name suggests, this example features a quadratic function given by ...
% f: R^2 -> R with 0 = f(z,mu) = mu - b + a * z^2.   (1)
% This equation is very similar to the saddle-node / fold bifurcation, but features the additional parameters a and b. The ...
% variable z corresponds to the state space vector. In this case, it is a scalar, but it is a vector in general. mu is called ...
% the continuation parameter and a and b are parameters influencing the shape and orientation of the parable. Here, let us set them to
a = 1;
b = 1;
%
% Equation (1) defines a curve in 2D space. Our aim is to plot the curve in a diagram where z is plotted against mu. CoSTAR will ...
% approximate discrete points on this curve by executing a path continuation. For a = b = 1, the plot will result in a parable ...
% which is opened to the left and which meets the horizontal axis at mu = 1.


%                 1.2 Continuation                  %
%
%          1.2.1 Setting useful variables           
%
% Before we specify the settings needed by CoSTAR, we introduce some important variables, which we will use later on. This is not ...
% necessarily needed but it helps to keep an overview of the most important variables/settings. Furthermore, important aspects ...
% of CoSTAR are explained here, making the following section "1.2.2 The CoSTAR settings" clearer.
IC = 0.1;                               % Sets the initial point z0 where the nonlinear systen solver fsolve starts to find the ...
%                                         first point on the curve. Should be as close as possible to the curve.
mu_limit = [-5, b];                     % Defines the limits of our continuation, i.e. the curve will be plotted within this interval.
mu0 = mu_limit(1);                      % Sets the mu-value at which the continuation begins.
%
Fcn = @(z,param) parable(z,param);      % Fcn is a function handle containing the right-hand side of equation (1).
% Fcn MUST ALWAYS be defined by "Fcn = @(z,param) ..." when calculating equilibrium solutions. param is a cell array storing all ...
% parameters, apart from z, that need to be passed to (1). The state variable vector z is always treated seperately and therefore ...
% not included in param. If (1) does not exhibit any parameters or all parameters are defined within Fcn ("parable" in this case), ...
% Fcn still needs to have the arguments (z,param) (this is due to the code being developed for continuation purposes primarily). ...
% Here, we need to define param since (1) exhibits the parameters a, b and mu.
%
% TASK: Open the function "parable" (click the right mouse button right on it and choose 'Open "parable"' or go into the "RHS" ...
%       subfolder located in the main path of CoSTAR and double click on the function). Have a look at how the equation is defined ...
%       and compare that to equation (1).
% IMPORTANT: Please pay attention to how the state variable z1 is defined and how mu, a and b are extracted from param.
%
% Now we define the "param" array. It is important that it is a cell array AND the order of the parameters placed in param must ...
% correspond to the definitions in the "parable" function! For example, the parameter a is defined by a = param{2} in the "parable" ...
% function. Therefore, a has to be the SECOND variable in the "param" array. For the variable mu, we use mu0 since mu0 defines the ...
% mu-value where the continuation starts.
param = {mu0, a, b};   
% Next, we have to tell CoSTAR where the continuation parameter is located within param. We also call this the "active" parameter. ...
% This is done by defining an integer variable which is used to get the continuation parameter from param by cell indexing, i.e. ...
% continuation parameter mu = param{active_parameter}. As we want to continue the curve via mu and the corresponding value mu0 is ...
% the FIRST element of param, we set:
active_parameter = 1;


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
% In this case, "options.opt_cont" is additionally needed.
% 
% Each of the "option" structures contain mandatory fields which ALWAYS have to be set. 
% Furthermore, there are optional fields that CAN be set. Most of the optional fields have a default value. 
% However, even optional fields NEED to be set in some cases (denoted as "Opt-Need.", e.g. 'param' in this example, see below).
% 
% IMPORTANT: All option structures must be created by the "costaropts" function. The syntax of the arguments of "costaropts" is ...
%            equivalent to the MATLAB struct function, i.e. costaropts(fieldname1,value1,...,fieldnameN,valueN).
%
% Let us start with the "options.system" structure:
options.system = costaropts('order',0,'dim',1,'rhs',Fcn,'param',param,'info','continuation of parable equation');
% Mandatory fields: - 'order':  Describes the order of the ODE. As we treat an equilibrium problem which is defined by the ...
%                               algebraic equation 0 = f(z, mu) (not an ODE), we set 'order' = 0.
%                   - 'dim':    Dimension of the state space. In this case, the state variable z is a scalar, so the dimension is 1.     
%                   - 'rhs':    Right-hand side of our equation, which we already assigned to the variable Fcn.           
% Opt-Need. field:  - 'param':  This field expects the previously defined "param" array. We NEED to set it because it is required ...
%                               by Fcn. Moreover, when executing a continuation, param is always necessary. However, param CAN be ...
%                               omitted when only calculating the initial solution (i.e. no continuation is carried out).
% Optional field:   - 'info':   Contains your individual text about the computation. This does not need to be set. However, it may ...
%                               help you later on if you have forgotten what this computation is about (default: no text).
%
% Next, there is the "options.opt_sol" structure:
options.opt_sol = costaropts('sol_type','equilibrium','cont','on','stability','off','act_param',active_parameter);
% Mandatory fields: - 'sol_type':   We need to specify what type of solution CoSTAR has to calculate. In this case, we want to ...
%                                   calculate an equilibrium solution. It is also possible to use 'eq' instead of 'equilibrium'.
%                   - 'cont':       As we want to plot a curve, we need to do a continuation, so we set 'cont' = 'on'. If ...
%                                   'cont' = 'off', CoSTAR would only calculate one single point.
%                   - 'stability':  CoSTAR computes the stability of the solution if 'stability' = 'on'. Stability is addressed in a ...
%                                   separate tutorial in detail. In order to keep this tutorial simple, we set the field to 'off'.
% Opt-Need. field:  - 'act_param':  As we want to do a continuation, we have to tell CoSTAR where the continuation parameter is ...
%                                   located within param. To do that, we use the previously definied variable "active_parameter". ... 
%                                   Here, this field NEEDS to be set, because we provided a "param" array. If there is no ...
%                                   "param" array, 'act_param' may not be set.
%
% Going on, we have to set the "options.opt_init" structure. The fields of the "options.opt_init" structure depend on the solution ...
% type as well as the chosen approximation method. For equilibrium solutions, there is only one mandatory field and there are ...
% no optional fields.
options.opt_init = costaropts('ic',IC);
% Mandatory fields: - 'ic': This property sets the initial point which is needed by the nonlinear systen solver fsolve in order ...
%                           to calculate the first point on the curve. Should be as close as possible to the curve. We already ...
%                           defined IC above.
%                           -> Allowed values: [dim x 1] (double) array             (no default value)
%
% So far we have definied all "options" structures which CoSTAR always needs.
% Since we want to do a continuation, we also have to set the "options.opt_cont" structure:
options.opt_cont = costaropts('mu_limit',mu_limit,'step_control','off');
% Mandatory fields: - 'mu_limit':       Sets the limits of the continuation. For this purpose, we defined the "mu_limit" variable.
% Optional fields:  - 'step_control':   Specifies the step control method to be used for the continuation. Using step control ...
%                                       is addressed in a separate tutorial. In order to keep this tutorial simple, we do not use ...
%                                       a step control and therefore set the field to 'off' (default: 'angle'). The step width will ...
%                                       be 0.1 consistently.
%
% Finally, we are done defining the required settings. All solution type specific fields, which are the fields of "options.opt_init", ...
% were explained above. Concerning the rest of the "options" stuctures ("options.system", "options.opt_sol" and "options.opt_cont"), ...
% only the necessary fields were defined and explained to some extend.
% If you want to have a deeper insight into the "options" structures and its fields, please use the "costarhelp" function by typing ...
% "costarhelp.options" in the command window. In order to directly open the help pages of particular "options" structures, type ...
% "costarhelp.<name_of_options_structure>", e.g. "costarhelp.opt_cont". In case of "options.opt_init", type ...
% "costarhelp.opt_init('EQ')" or "costarhelp.opt_init('equilibrium')".


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
% Apart from the figure plot, CoSTAR displays information in the command window during the computation:
% - At the beginning, the iteration process of fsolve trying to find the first point on the curve is shown.
% - As soon as fsolve succeded, CoSTAR announces "Initial solution found!".
% - After that, CoSTAR displays "Iter: <XXX> -- mu = <XXX> -- stepwidth = <XXX>" when a new point on the curve has been calculated.
%       * "Iter" depicts the number of points on the curve which already have been calculated.
%       * "mu" shows the value of the continuation parameter mu at the latest point on the curve.
%       * "stepwidth" displays the step width that was used to calculate the latest point on the curve.
% - Finally, CoSTAR reports the reason of termination of the continuation.


%               1.3  Postprocessing                 %
%
% As you might have noticed, CoSTAR plots the maximum of the Euclidean norm of the state space vector (which is equivalent to ...
% abs(z) in this case). That is why we do not see the parable yet. In order to change the vertical axis and to plot z against mu, ...
% we need to call the "contplot" function. Similar to the "costar" function, contplot expects a structure that defines all required ...
% options for the plot. 
opt_contplot = costaropts('zaxis', @(z) max(z(:,1)));
% By setting the mandatory field 'zaxis' to "@(z) max(z(:,1))", only the maximum of the first element of the state space vector ...
% (which is equivalent to z itself in this case) is plotted on the vertical axis. You might ask yourself why we need the "max" ...
% function and can not define costaropts('zaxis', @(z) z(:,1)). The reason for this is that the "contplot" function (see below) ...
% also has to work when computing periodic and quasi-periodic solutions, where "max(z(:,1))" is in fact required.
%
% Now we can call the "contplot" function which creates a new continuation plot. It is a function of the solution class object and ...
% therefore has to be called by "S.contplot(...)". Apart from the "opt_contplot" structure, contplot also requires the ...
% DynamicalSystem class object DYN. contplot returns the calculated z and mu values.
[z,mu] = S.contplot(DYN,opt_contplot);              % Create a new continuation plot
%
% Further postprocessing options and functions can be found via the costarhelp function (type "costarhelp.costar" in the command ...
% window) since a complete explanation of these would exceed the scope of this tutorial.


% There's the parable and we are done with our first example! Now we move on to the second example "pitchfork bifurcation".



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Ex. No 2: Pitchfork Bifurcation        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear variables; clc; close all;                    % Let's clean our "desk" before we begin with our second example.


%                2.1 Introduction                   %
%
% The second example we have a look on is the pitchfork bifurcation which is given by
% f: R^2 -> R with 0 = f(z,mu) = mu * z - z^3.   (2)
% Again, z corresponds to the state variable and mu is the continuation parameter. The function f is given in a very similar way ...
% in CoSTAR, namely by f(z,mu) = mu * z - z^3 + gamma. Obviously, f exhibits the additional parameter gamma that needs to be defined ...
% by the user. For gamma = 0, the function f describes the pitchfork bifurcation, having a bifurcation point at mu = z = 0.
gamma = 0;
% (The implemented function exhibits gamma so that it can further be used to calculate imperfections of the pitchfork bifurcation.)
%
% As in the "Parable" example, we want to compute the solution set of (2) by executing a path continuation and plot z against mu.
% This time, however, you will be asked to apply some of the things that you have just learnt about CoSTAR!
%
% NOTE: In the following it is assumed that you have gone through the "Parable" example as well as that you know how the ...
%       pitchfork bifurcation looks like! If you have not worked through the "Parable" example yet, please do so before ...
%       continuing with this example. If you do not know how the pitchfork bifurcation looks like, place the cursor in this ...
%       section of the script and click "Run Section". This will generate a plot in which the pitchfork bifurcation is shown.


%             2.2 Continuation - Part 1             %
%
%          2.2.1 Setting useful variables           
%
mu_limit = [-1.5, 1.5];                             % This time, we set the continuation limits to [-1.5, 1.5].
mu0 = mu_limit(1);                                  % Like before, we want to begin the continuation at the bottom limit of mu.
IC = 0;                                             % We already know that z = 0 is a solution of f(z,mu) = 0, so we set the ...
%                                                     initial value to 0. As a result, fsolve does not need to iterate to the ...
%                                                     first point on the curve.                                                  
Fcn =  @(z,param) pitchfork_ap(z,param);            % Fcn contains the right-hand side of 0 = f(z,mu).
%
% TASK: Open the function "pitchfork_ap" and have a look at the definitions. 
% QUESTION: How does the "param" array need to be defined and what is the value of the variable "active_parameter"? 
%           (The correct answers can be found at the end of the following lines.)
% param = ???                                                                                                                                                               
% active_parameter = ???                                                                                                                                                    
                                                                                                                                                        param = {mu0, gamma};
                                                                                                                                                        active_parameter = 1;
%            2.2.2 The CoSTAR settings              
%
% Now we can set the "options" structures.
% QUESTIONS: - Which "options" structures does CoSTAR require in any case? (Answers to the right -> )                                                   "options.system", "options.opt_sol", "options.opt_init"
%            - Which "options" structures do we need additionally?                                                                                                                           "options.opt_cont"
%
options1.system   = costaropts('order',0,'dim',1,'rhs',Fcn,'param',param,'info','continuation of pitchfork bifurcation - part 1');
options1.opt_sol  = costaropts('sol_type','equilibrium','cont','on','stability','on','act_param',active_parameter);
options1.opt_init = costaropts('ic',IC);
options1.opt_cont = costaropts('mu_limit',mu_limit,'step_control','off');
% Overall, the settings are fairly identical to the setting of the "Parable" example.
% We know that one branch will change its stability behaviour during the continuation. In order to detect that, we set the field ...
% 'stability' of the "options.opt_sol" structure to 'on'. (For more information on calculating the stability of solutions, please see ...
% the Stability Tutorial.) Again, if you want to have a deeper insight into the "options" structures and their fields, please use ...
% the "costarhelp" function.


%     2.2.3 Calling CoSTAR and running the simulation
%
% In the next step, we can already call CoSTAR.
[S1,DYN1] = costar(options1);                       % Calling CoSTAR and performing the continuation
%
% As you can see, CoSTAR is not able to automatically continue multiple branches yet when running into a bifurcation point. Only ...
% the solution z = 0 is displayed. For mu < 0, the blue line indicates that the solution is stable, while for mu > 0, the red line ...
% shows that the solution is unstable. At mu = 0, CoSTAR has detected a bifurcation point, which is shown in the plot ("FB": fold, ...
% pitchfork or transcritical bifurcation). Details are saved to the field 'bifurcation' of the solution object S1.
% It is pointed out again that details of the computation of stability and bifurcation points are addressed in a separate tutorial ...
% and are therefore not explained here.
%
% So far, CoSTAR has calculated the solution z = 0, but we also want to compute the top and bottom branches for mu >= 0. For this ...
% purpose, we need to do a second path continuation.


%             2.3 Continuation - Part 2             %
%
% Most of the settings needed for the second continuation can be taken from the continuation above. 
% However, we need to make some adjustments in order to compute the missing branches via a single continuation.
%
% QUESTION: How does mu0 and IC need to be defined in order to compute the missing branches via a single continuation?
%           (Possible correct answers can be found at the end of the following lines.)
% mu0 = ???
% IC = ???
                                                                                                                                                        mu0 = mu_limit(2);
                                                                                                                                                                   IC = 1;
% Explanation of mu0: As we want to compute the missing branches via a single continuation, we obviously need to start at ...
%                     the upper limit of the continuation parameter, which is why we set mu0 = mu_limit(2).
% Explanation of IC: The initial point for fsolve to find the first point on the branch must be close enough to the branch ...
%                    so that fsolve converges to the desired branch. For example, when setting IC = 0.5, fsolve converges ...
%                    to the solution branch at z = 0.
%
% Apart from mu0, IC and the "info" field, we need to make two more adjustments to the "options". First, we set the continuation ...
% step width to 0.05 (default: 0.1). This is not required but it ensures a good resolution of the curve. Second, we need to change ...
% the direction of the continuation at the start by setting the field 'direction' of the "opt_cont" structure to -1. This is ...
% required here, because CoSTAR always starts to continue the curve in positive mu direction by default ('direction' = 1). However, ...
% CoSTAR must begin this continuation in negative mu direction since we start at the upper limit of mu. You can try to not set ...
% 'direction' to -1, but the continuation would stop immediately (pay attention to the termination message).
%
param = {mu0, gamma};                               % Since mu0 changed, we need to update param.
options2.system   = costaropts('order',0,'dim',1,'rhs',Fcn,'param',param,'info','continuation of pitchfork bifurcation - part 2');  
options2.opt_sol  = costaropts('sol_type','equilibrium','cont','on','stability','on','act_param',active_parameter);
options2.opt_init = costaropts('ic',IC);
options2.opt_cont = costaropts('mu_limit',mu_limit,'step_control','off','step_width',0.05,'direction',-1);
%
% Now we have set the required options to calculate the missing branches. 
[S2,DYN2] = costar(options2);                       % Calling CoSTAR and performing the second continuation


%                2.4 Postprocessing                 %
%
% All calculations are done now and we see the plots of our two continuations. Similar to the "Parable" example, we need to create ...
% a new plot in order to depict z (and not norm(z)) on the vertical axis. We already know how to do that and call the "contplot" ...
% function of the solution object S1. This will create a plot in which the solution branch z = 0 is displayed in a diagram ...
% depicting z on the vertical axis.
opt_contplot1 = costaropts('zaxis', @(z) max(z(:,1)));
[z1,mu1] = S1.contplot(DYN1,opt_contplot1);
%
% In contrast to the "Parable" example, we need to complement the plot by the second continuation in order to depict all branches ...
% within one diagram. This is done by adding the field 'figure' to the "opt_contplot" structure and setting 'figure' to the current ...
% figure handle (gcf). As the curve of the second continuation is missing, we call the "contplot" function of S2.
opt_contplot2 = costaropts('zaxis', @(z) max(z(:,1)), 'figure', gcf);
[z2,mu2] = S2.contplot(DYN2,opt_contplot2);
%
% Now all solution branches are shown in one diagram. Unfortunately, the last contplot call displaced the horizontal axis limits. ...
% Therefore, we rescale the limits of the horizontal axis to the limits of the continuation parameter mu.
xlim(mu_limit)
%
% Further postprocessing options and functions can be found via the costarhelp function (type "costarhelp.costar" in the command ...
% window) since a complete explanation of these would exceed the scope of this tutorial.

% Done!



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Final Words                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The CoSTAR tutorial on calculating equilibrium solutions is now finished.
% For additional information, please use the "costarhelp" function and/or the CoSTAR manual.

% If you are interested in learning about further capabilities of CoSTAR, you are invited to have a look at the other tutorials as well.

% See you soon!