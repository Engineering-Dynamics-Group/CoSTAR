% Gatekeeper function for the SHM subclass. In here, all input parameters are checked, before further processing. 
% This static method is called from the static method AM_gatekeeper.
%
% @GC:                  object of Gatekeeper class
% @system:              user supplied option structure for the system
% @opt_approx_method:   user supplied option structure for the solution method
% @opt_init:            user supplied option structure for the initial condition

function s_PS_SHM_gatekeeper(GC,system,opt_sol,opt_approx_method,opt_init)


%% Check the opt_approx_method structure
opt_approx_mandatory_fieldnames  = {};                                  % mandatory fieldnames in the options structure
opt_approx_allowed_fieldnames    = {'solver','n_shoot','phase_condition'};   % allowed fieldnames in the options structure

GC.check_fields(opt_approx_method,'opt_approx_method',opt_approx_mandatory_fieldnames,opt_approx_allowed_fieldnames);
GC.speak;


%% Check the opt_approx_method field entries

% Check the mandatory fields first (these are definitively present)
%%%%%%%%%%%%%%%%%%%%

% Check the optional fields now
%%%%%%%%%%%%%%%%%%%%
solver_allowed_fieldvalues = {'ode45','ode78','ode89','ode23','ode113','ode15s','ode23s','ode23t','ode23tb'};
if isfield(opt_approx_method,'solver')
    GC.check_data(opt_approx_method.solver,'opt_approx_method.solver','char',[],solver_allowed_fieldvalues); 
end
GC.speak;

if isfield(opt_approx_method,'n_shoot')            
    GC.check_data(opt_approx_method.n_shoot,'opt_approx_method.n_shoot','double','scalar','positive');
    if mod(opt_approx_method.n_shoot,1)~=0
        GC.error_msg{1,end+1} = append('The value of the options field opt_approx_method.n_shoot is ',num2str(opt_approx_method.n_shoot),', but it must be an integer!');
    end
    if opt_approx_method.n_shoot == 0
        GC.error_msg{1,end+1} = append('The value of the options field opt_approx_method.n_shoot is 0, but it must be a positive integer (> 0)!');
    end
end
GC.speak;

phase_condition_allowed_fieldvalues = {'poincare','integral'};
if isfield(opt_approx_method,'phase_condition')
    GC.check_data(opt_approx_method.phase_condition,'opt_approx_method.phase_condition','char',[],phase_condition_allowed_fieldvalues);
    if isfield(opt_sol,'non_auto_freq')
        GC.error_msg{1,end+1} = 'You defined a phase condition via opt_approx_method.phase_condition and';
        GC.error_msg{1,end+1} = 'the system seems to be non-autonomous since you set the option opt_sol.non_auto_freq.';
        GC.error_msg{1,end+1} = 'However, a phase condition is only allowed for autonomous systems.';   
    end
end



%% Check the opt_init structure
opt_init_mandatory_fieldnames  = {'ic'};                                % mandatory fieldnames in the options structure
opt_init_allowed_fieldnames    = {'ic'};                                % allowed fieldnames in the options structure

GC.check_fields(opt_init,'opt_init',opt_init_mandatory_fieldnames,opt_init_allowed_fieldnames);
GC.speak;


%% Check the opt_init field entries

%Check the mandatory fields first (these are definitively present)
%%%%%%%%%%%%%%%%%%%%
GC.check_data(opt_init.ic,'opt_init.ic','double','vector',[]);
GC.speak;
if ~isfield(opt_approx_method,'n_shoot');       n_shoot = 5;
else;                                           n_shoot = opt_approx_method.n_shoot;
end
if ~((numel(opt_init.ic) == system.dim) || (numel(opt_init.ic) == system.dim*n_shoot))
    GC.error_msg{1,end+1} = append('The options field opt_init.ic is a ',num2str(numel(opt_init.ic)),'-dimensional vector.');
    GC.error_msg{1,end+1} = 'However, the length of opt_init.ic must be equal either to the dimension dim of the system or to dim*n_shoot, ';
    GC.error_msg{1,end+1} = 'where n_shoot is the number of shooting points.';
end
GC.speak;


% Check the optional fields now
%%%%%%%%%%%%%%%%%%%%


end