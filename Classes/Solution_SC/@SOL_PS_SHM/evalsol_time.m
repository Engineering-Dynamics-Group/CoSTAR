% Method of SOL_PS_MSHM: This method (re-)calculates a periodic orbit of a multiple shooting solution
% starting from the iterated point on the periodic oribt
%
% @obj:     Solution subclass object
% @DYN:     DynamicalSystem object
% @options: options structure for postprocessing solutions
%
% @s:       Time solution array: This must(!) be a [options.resolution x state_space_dimension x n_evals] dimensional array !!!
% @mu:      Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @t:       Array of the time points: This must(!) be a [options.resolution x 1 x n_evals]  dimensional array !!!
% n_evals:  Number of curve points to be evaluated 

function  [s,mu,t] = evalsol_time(obj,DYN,options)

    index = options.index;
    dim = DYN.dim;
    N = numel(index);
  
    for k =1:N

        idx = index(k);

        freq    = obj.freq(1,idx);
        mu      =  obj.mu(1,idx);
        s0      =  obj.s(:,idx);

        if isfield(options,'interval')          %If an integration interval was supplied by the user... use this interval - else: used the standard interval
            tspan = linspace(options.interval(1),options.interval(2),options.resolution);
        else
            tspan = linspace(0,2.*pi./freq,options.resolution);                                                     %Integration interval
        end

        Fcn = DYN.rhs;
        param = DYN.param;                                                                                  %Set parameters
        param{DYN.act_param} = mu;
        
        if tspan(1) > 0                                                                                     %Start-up integration necessary ...
            [~,s_su] = obj.solver_function(@(t,z)Fcn(t,z,param),[0,tspan(1)],s0(1:dim,1),obj.odeOpts);      %in the range of [0,tspan(1)]
            s1 = s_su(end,:)';                                                                              %Starting point of actual integration
        else                                                                                        
            s1 = s0;                                                                                        %No start-up integration necessary
        end

        [t(:,1,k) ,s(:,:,k) ] = obj.solver_function(@(t,z)Fcn(t,z,param),tspan,s1(1:dim,1),obj.odeOpts);    % Start evaluate from first shooting point and integrate over whole interval
        
    end

    %Get the mu values
    mu = obj.mu(options.index);    %options.index is unique due to S.solget_up_index

end