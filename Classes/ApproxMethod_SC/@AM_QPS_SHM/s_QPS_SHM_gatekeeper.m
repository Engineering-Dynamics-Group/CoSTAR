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
opt_init_allowed_fieldnames    = {'iv','c0','c1_matrix','s1_matrix'};       %allowed fieldsnames in the opt_init structure

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
if((isfield(opt_init,'c0')||isfield(opt_init,'c1')||isfield(opt_init,'s1'))&&isfield(opt_init,'iv'))
    GC.error_msg{1,end+1} = 'You supplied a fourier series and "iv", but you must supply only one or none!';
end

% If no initial conditions are provided, set fourier-series to zeros
if(~(isfield(opt_init,'c0')||isfield(opt_init,'c1')||isfield(opt_init,'s1')||isfield(opt_init,'iv')))
    opt_init.c0 = zeros(system.dim,1); 
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

% Field 'c0'
if isfield(opt_init,'c0')
    GC.check_data(opt_init.c0, 'opt_init.c0', 'double', [system.dim,1], []);
    GC.speak();
end

% Field 'c1_matrix'
if isfield(opt_init,'c1_matrix')
    GC.check_data(opt_init.c1_matrix, 'opt_init.c1_matrix', 'double', [], []);
    if size(opt_init.c1_matrix,1) ~= system.dim         % Check that the number of rows matches the system dimension
        GC.error_msg{1,end+1} = append('The number of rows of the opt_init structure field "c1_matrix" is ', num2str(size(opt_init.c1_matrix,1)), '.');
        GC.error_msg{1,end+1} = append('However, the number of rows must match the dimension of the system, which is ', num2str(system.dim), '.');
    end
    if size(opt_init.c1_matrix,2) > 3                   % The number of columns can be <= 3 (missing columns are filled by zero vectors)
        GC.error_msg{1,end+1} = append('The number of columns of the opt_init structure field "c1_matrix" is ', num2str(size(opt_init.c1_matrix,2)), '.');
        GC.error_msg{1,end+1} = 'However, the number of columns must not be greater than 3.';
    end
    GC.speak();
end

% Field 's1_matrix'
if isfield(opt_init,'s1_matrix')
    GC.check_data(opt_init.s1_matrix, 'opt_init.s1_matrix', 'double', [], []);
    if size(opt_init.s1_matrix,1) ~= system.dim         % Check that the number of rows matches the system dimension
        GC.error_msg{1,end+1} = append('The number of rows of the opt_init structure field "s1_matrix" is ', num2str(size(opt_init.s1_matrix,1)), '.');
        GC.error_msg{1,end+1} = append('However, the number of rows must match the dimension of the system, which is ', num2str(system.dim), '.');
    end
    if size(opt_init.s1_matrix,2) > 3                   % The number of columns can be <= 3 (missing columns are filled by zero vectors)
        GC.error_msg{1,end+1} = append('The number of columns of the opt_init structure field "s1_matrix" is ', num2str(size(opt_init.s1_matrix,2)), '.');
        GC.error_msg{1,end+1} = 'However, the number of columns must not be greater than 3.';
    end
    GC.speak();
end

if isfield(opt_approx_method,'solver'); GC.check_data(opt_approx_method.solver,'opt_approx_method.solver','char',[],solver_allowed_fieldvalues); end
if isfield(opt_approx_method,'n_char'); GC.check_data(opt_approx_method.n_char,'opt_approx_method.n_char','double','scalar',[]); end
GC.speak;

end