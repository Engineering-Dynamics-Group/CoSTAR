% Gatekeeper function for the SHM subclass. In here, all input parameters are checked, before further processing. 
% This static method is called from the static method AM_gatekeeper.
%
% @GC:                  object of Gatekeeper class
% @system:              user supplied option structure for the system
% @opt_approx_method:   user supplied option structure for the solution method
% @opt_init:            user supplied option structure for the initial condition

function s_PS_SHM_gatekeeper(GC,system,opt_approx_method,opt_init)


%% Check the opt_approx_method structure
opt_approx_mandatory_fieldnames  = {};                                  % mandatory fieldnames in the options structure
opt_approx_allowed_fieldnames    = {'solver','n_shoot'};                % allowed fieldnames in the options structure

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



%% Check the opt_init structure
opt_init_mandatory_fieldnames  = {'ic'};                                % mandatory fieldnames in the options structure
opt_init_allowed_fieldnames    = {'ic'};                                % allowed fieldnames in the options structure

GC.check_fields(opt_init,'opt_init',opt_init_mandatory_fieldnames,opt_init_allowed_fieldnames);
GC.speak;


%% Check the opt_init field entries

%Check the mandatory fields first (these are definitively present)
%%%%%%%%%%%%%%%%%%%%
GC.check_data(opt_init.ic,'opt_init.ic','double',[system.dim,1],[]);
GC.speak;

% Check the optional fields now
%%%%%%%%%%%%%%%%%%%%


end