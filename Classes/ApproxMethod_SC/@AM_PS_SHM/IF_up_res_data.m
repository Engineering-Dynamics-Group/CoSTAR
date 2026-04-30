% This function is a method of the subclass AM_PS_SHM.
% It passes information between the continuation algrithm and the ApproxMethod subclass.
%
% @obj: ApproxMethod subclass AM_PS_SHM object
% @CON: Continuation class object

function obj = IF_up_res_data(obj,var1,DYN)

    if DYN.n_auto == 1                      % If the system is autonomous        
        
        %% Get the correct initial value
        if isa(var1,'Continuation')                 % If var1 is an object of Continuation
    
            switch obj.phase_condition              % Update the current initial condition. Used for the phase condition
                case 'poincare'
                    obj.iv = var1.yp(1:end-1);      % Take the predictor point
                    mu = var1.yp(end);              % Get the corresponding mu-value
                case 'integral'
                    obj.iv = var1.y0(1:end-1);      % Take the preceding solution
                    mu = var1.y0(end);              % Get the corresponding mu-value
            end
    
        elseif isa(var1,'double')                   % var1 should be a solution vector (type double) in all other cases
    
            obj.iv = var1(1:end-1);                 % Set iv to given solution vector x0 (only relevant in initial_solution at the error control)
            mu = var1(end);                         % Get the corresponding mu-value
    
        end   
        
        
        %% If the integral phase condition is used: Integrate the preceding solution to evaluate the phase condition
        if strcmpi(obj.phase_condition,'integral')

            % Parameters
            dim = DYN.dim;                      % Dimension of the system
            Fcn = DYN.rhs;                      % RHS of ODE
            param = DYN.param;                  % Parameter array
            n_shoot = obj.n_shoot;              % Number of shooting points
            n_time = obj.n_time;                % Number of time evaluation points in each shooting interval for the integral phase condition
            
            x0_mat = reshape(obj.iv(1:end-1),dim,n_shoot);              % Get the shooting points of iv and reshape them to a matrix
            dT0 = 2*pi/(obj.iv(end)*n_shoot);                           % Time interval between two shooting points (the last element of iv is omega)
            T0_int = [0:dT0:(n_shoot-1)*dT0; dT0:dT0:n_shoot*dT0].';    % Start and end times for the integration
            
            % calc_stability = strcmpi(DYN.stability,'on');               % Logical variable defining if stability is computed
            % if calc_stability;  h = eps^(1/3);                          % When stability is computed: Jacobian is approximated using central finite difference
            % else;               h = sqrt(eps);                          % When stability is not computed: Jacobian is approximated using forward finite difference
            % end
            
            % Initialisations
            Z0 = NaN(dim,n_time,n_shoot);                               % Initialisation
            Z0_mu_plus = Z0;                                            % Initialisation
            Z0_mu_minus = Z0;                                           % Initialisation
            param{DYN.act_param}  = mu;                                 % Update the active parameter
            % param_mu_plus = param;                                      % Define a new parameter array
            % param_mu_minus = param;                                     % Define a new parameter array
            % param_mu_plus{DYN.act_param}  = mu + h*(1+abs(mu));         % Update the new parameter array by mu + delta_mu
            % param_mu_minus{DYN.act_param} = mu - h*(1+abs(mu));         % Update the new parameter array by mu - delta_mu
            
            % Integration
            for k = 1:n_shoot
                [~,Z] = obj.solver_function(@(t,z) Fcn(t,z,param), linspace(T0_int(k,1),T0_int(k,2),n_time+1), x0_mat(:,k), obj.odeOpts);
                Z0(:,:,k) = Z(1:end-1,:).';                             % Save the trajectory (not Z(t_end) because t_end already belongs to the next interval)
                % Now do a second/third integration with perturbed mu-value (needed for Jacobian)
                % if ~(isfield(DYN.system,'first_integral') && (DYN.act_param == numel(param)))  % If the value of the first integral is NOT the continuation parameter (if it is: dg/dmu = 0)
                %     [~,Z_mu_plus]  = obj.solver_function(@(t,z) Fcn(t,z,param_mu_plus), linspace(T0_int(k,1),T0_int(k,2),n_time+1), x0_mat(:,k), obj.odeOpts);
                %     Z0_mu_plus(:,:,k) = Z_mu_plus(1:end-1,:).';         % Save the "+ delta" perturbed mu trajectory at the predictor point
                %     if calc_stability                                   % Only required when using central finite difference for Jacobian
                %         [~,Z_mu_minus] = obj.solver_function(@(t,z) Fcn(t,z,param_mu_minus), linspace(T0_int(k,1),T0_int(k,2),n_time+1), x0_mat(:,k), obj.odeOpts);
                %         Z0_mu_minus(:,:,k) = Z_mu_minus(1:end-1,:).';   % Save the "- delta" perturbed mu trajectory at the predictor point
                %     end
                % end
            end
            
            % Save the trajectory in obj
            obj.Z0{1} = Z0;                                             
            obj.Z0{2} = Z0_mu_plus;
            obj.Z0{3} = Z0_mu_minus;
            
        end

    end

end