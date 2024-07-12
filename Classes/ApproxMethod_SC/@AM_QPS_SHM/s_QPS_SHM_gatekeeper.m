%Gatekeeper function for the AM_QPS_SHM subclass. In here, all input
%parameters are checked, before further processing. This static method is
%called from the static method AM_gatekeeper
%
%@GC:                   object of Gatekeeper class
%@system                user supplied option structure for the DynamicalSystem class                   
%@opt_sol_method:       user supplied option structure for the solution method                  
%@opt_approx_method:    user supplied option structure for the approximation method                   
%@opt_init:             user supplied option structure for the initial condition
                       
function s_QPS_SHM_gatekeeper(GC,system,opt_sol,opt_approx_method,opt_init)

opt_approx_method_mandatory_fieldnames  = {};                               %mandatory fieldsnames in the opt_approx structure
opt_approx_method_allowed_fieldnames    = {'solver','n_char'};              %allowed fieldsnames in the opt_approx structure

opt_init_mandatory_fieldnames  = {};                                        %mandatory fieldsnames in the opt_init structure
opt_init_allowed_fieldnames    = {'iv','ic','tinit','deltat','dt'};         %allowed fieldsnames in the opt_init structure

%% Check the opt_approx_method ans opt_init_method structure
GC.check_fields(opt_approx_method,'opt_approx_method',opt_approx_method_mandatory_fieldnames,opt_approx_method_allowed_fieldnames);
GC.check_fields(opt_init,'opt_init',opt_init_mandatory_fieldnames,opt_init_allowed_fieldnames);
GC.speak();

%% Check the entries
solver_allowed_fieldvalues = {'ode45','ode78','ode89','ode23','ode113','ode15s','ode23s','ode23t','ode23tb'};
%Check the mandatory fields first (these are defintively present)
%%%%%%%%%%%%%%%%%%%%

%Check the optional fields now
%%%%%%%%%%%%%%%%%%%%
if(isfield(opt_init,'ic')&&isfield(opt_init,'iv'))
    GC.error_msg{1,end+1} = 'You supplied "ic" and "iv", but you must supply only one or none!';
end

if(isfield(opt_init,'ic'))
    if(~isempty(opt_init.ic))
        GC.check_data(opt_init.ic,'opt_init.ic','double','vector',[]);
    else
        opt_init.ic = zeros(system.dim,1);
    end
    if(size(opt_init.ic,1)~=system.dim)
        GC.error_msg{1,end+1} = append('You supplied "ic" which has size ',num2str(size(opt_init.ic,1)),' x 1 but should be ',num2str(system.dim),' x 1!');
    end
end

if(isfield(opt_init,'iv'))
    GC.check_data(opt_init.iv,'opt_init.iv','double','vector',[]);
    if mod(numel(opt_init.iv)/system.dim,1) ~= 0        % n_char_iv = numel(opt_init.iv)/system.dim must be an integer
        GC.error_msg{1,end+1} = append('You supplied "iv", which has size [',num2str(size(opt_init.iv,1)),' x 1]. However, the number of characteristics of "iv" equals');
        GC.error_msg{1,end+1} = append('numel(iv)/system.dim = ', num2str(numel(opt_init.iv)/system.dim), ', which is not an integer. This is not allowed.');
        GC.error_msg{1,end+1} = 'If this error occurs, please make sure that numel(iv) is a product of system.dim and another integer.';
    end
    % Not needed anymore:
    % if(size(opt_init.iv,1)~=(system.dim*opt_approx_method.n_char))
    %     GC.error_msg{1,end+1} = append('You supplied "iv" which has size ',num2str(size(opt_init.iv,1)),' x 1 but should be ',num2str(system.dim*opt_approx_method.n_char),' x 1!');
    % end
end

if(isfield(opt_init,'tinit'))
    GC.check_data(opt_init.tinit,'opt_init.tinit','double','scalar',[]);
end
if(isfield(opt_init,'deltat'))
    GC.check_data(opt_init.deltat,'opt_init.deltat','double','scalar',[]);
end
if(isfield(opt_init,'dt'))
    GC.check_data(opt_init.dt,'opt_init.dt','double','scalar',[]);
end

if(isfield(opt_init,'deltat')&&isfield(opt_init,'dt'))                      % Check if dt is smaller than the integration time
    if(opt_init.dt>opt_init.deltat) 
        GC.error_msg{1,end+1} = append('The integration time on the attractor deltat = ',num2str(opt_init.deltat),' is shorter than your increment dt = ',num2str(opt_init.dt),'. Choose a lower dt!');
    end
end

if isfield(opt_approx_method,'solver'); GC.check_data(opt_approx_method.solver,'opt_approx_method.solver','char',[],solver_allowed_fieldvalues); end
if isfield(opt_approx_method,'n_char'); GC.check_data(opt_approx_method.n_char,'opt_approx_method.n_char','double','scalar',[]); end
GC.speak;


end






