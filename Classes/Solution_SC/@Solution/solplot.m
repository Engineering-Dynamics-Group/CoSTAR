% This function plots prescribed solution points in a 2D or 3D plot
%
% Input Arguments
% @obj:                 Object of Solution class
% @DYN:                 Object of DynamicalSystem class
% @options:             struct contain options for calculating solution
%
%                       Mandatory fields:
%                       'zaxis':    -->'function_handle:' User supplied
%                                       function handle. Example: options('zaxis',@(z)z(:,1))
%                                       for options.space = 'trajectory'
%                                       this is the only option.
%                                   --> 'euclidean': gives back the
%                                       euclidean (L2) along the states of an output vector
%                                       for the specified index.
%                                       Example: periodic solution
%                                       options('zaxis','euclidean') for an
%                                       200 x 3 x 10 (200 points per limit cycle, 3 states, 10 conti points)
%                                       a 200 x 1 x 10 vector
%                                   --> 'all': is not altering the output
%
%                       'space':    --> 'trajectory': requires 'zaxis' and
%                                       'xaxis'. Both must be function
%                                       handles. 'yaxis' for 3D plot
%                                       optional. Plot the values of the
%                                       function handles over each other
%                                   --> 'time': plots 'zaxis' over time.
%                                   --> 'hypertime: plot 'zaxis' over one or
%                                       two hypertime axis (2D/3D) Periodic:
%                                       2D time trajectory over torus
%                                       coordinate \theta = [0,2\pi]
%                                       Quasi-periodic: 3D plane over
%                                       torus coordinate \theta_1,2 =
%                                       [0,2\pi]x [0, 2\pi]
%                                   --> 'frequency': plots the FFT of
%                                       the hypertime solution (1D or nD FFT)
%                                       over the Angular frequency
%
%                       Optional fields:
%                       'xaxis'         function handle. Only needed for
%                                       options.space = 'trajectory'
%
%                       'zaxis'         function handle. Only needed for
%                                       options.space = 'trajectory
%
%                       'figure'        handle to a existing figure. If
%                                       supplied, solplot plots into the
%                                       existing figure environment.
%
%                       'resolution'    1x1 or 1xnumb_base_frq double array
%                                       defining the number of points in
%                                       time or hypertime at which the
%                                       solution is evaluated. This is the
%                                       total number of points (see 'interval')
%                                       Default value: 200
%                                       Ignored for options.space = 'solution'
%
%                       'interval'      1x2 interval: Defines the start and
%                                       end point of the evaluation time
%                                       interval
%                                       Ignored for options.space =
%                                       'solution', 'hypertime'
%                                       Default value: [0,2\pi/frq]
%
%                       'index'         1xn array indicating the indices
%                                       of the solution curve, where the
%                                       solution is evaluated.
%                                       Default: All solutions are given
%                                       back
%                       'mu'            1xn array indicating the mu values
%                                       of the solution curve, where the
%                                       solution is evaluated. If the exact
%                                       value is not present, the closest
%                                       value is used. If these leads to
%                                       doubled values, they are ignored.
%                                       Mu might lead to ambiguities in
%                                       case of overhanging curves
%                       'color'         Define a color for the plot either by 
%                                       providing an rgb array or a char 
%                                       'r','g','b','c','m','y' or 'k'
%
%
% Output Arguments:
% @output:              Struct containing plotted solution data, domain data (if possible) and mu-values. 
%                       Fields depend on solution space:
%                        - options.space = 'time'
%                           output.z:    zaxis values
%                           output.time: time values
%                           output.mu:   mu values
%                        - options.space = 'trajectory' (2D)
%                           output.z:  zaxis values
%                           output.x:  xaxis values
%                           output.mu: mu values
%                        - options.space = 'trajectory' (3D)
%                           output.z:  zaxis values
%                           output.x:  xaxis values                           
%                           output.y:  yaxis values
%                           output.mu: mu values
%                        - options.space = 'hypertime'
%                           output.z:         zaxis values
%                           output.hypertime: hypertime values
%                           output.mu:        mu values
%                        - options.space = 'frequency'
%                           output.amplitude: absolute amplitude values
%                           output.angle:     phase angle values
%                           output.frequency: angular frequency values
%                           output.mu:        mu values
%
%
% Example:   options = costaropts('xaxis',@(z)z(:,1),'zaxis',@(z)z(:,2),'space','trajectory','mu',[0.5,1,1.5]);
%            output = S.solplot(DYN,options);


function output = solplot(obj,DYN,options)

    %% Gatekeeper check of input: Only the input not checked by solget_gatekeeper is checked here
    options = obj.solplot_gatekeeper(DYN,options);
    solget_options = options;                       % solget_options is a subset of solplot options.

    % Remove the fields from solget_options = (solplot) options, which are not allowed in solget - otherwise solget_gatekeeper complains
                                            solget_options = rmfield(solget_options,'zaxis');
    if isfield(solget_options,'xaxis');     solget_options = rmfield(solget_options,'xaxis');       end
    if isfield(solget_options,'yaxis');     solget_options = rmfield(solget_options,'yaxis');       end
    if isfield(solget_options,'figure');    solget_options = rmfield(solget_options,'figure');      end
    if isfield(solget_options,'color');     solget_options = rmfield(solget_options,'color');       end
    if isfield(solget_options,'linestyle'); solget_options = rmfield(solget_options,'linestyle');   end

    % Define the linestyle
    if isfield(options,'linestyle'); linestyle = options.linestyle; else; linestyle = '-'; end


    
    %% Plot

    %Why do we need this: originally fig = 1 was set, if nothing  was supplied.
    %However, calling solplot multiple times without specifying options.figure
    %led to plotting always to the same figure with the handle "1". You could
    %detect the current handle of all open figures and assign a different one to
    %fig, but that seemed more hacky than this one.

    if isfield(options,'figure')
        figure(options.figure);                     % Use existing figure
        set(gcf,'DefaultLineLineWidth',2)           
    else
        figure;
        set(gcf,'Color','w','NumberTitle','off','DefaultLineLineWidth',2,'Name','Plot by SOLUTION.solplot');
        set(gca,'TickLabelInterpreter','latex');    % Set tick label interpreter to LaTeX
        grid on;
    end

    hold on;                                                    % Allow multiple lines to be plotted and enable grid
    counter = 0;                                                % Counter needed for legend
    LegStr_old = get(findobj(gcf,'type','Legend'),'String');    % Get existing legend. Output is empty or a [1 x n] cell (n: number of legend entries)
    sol_type = char(DYN.sol_type);                              % DYN.sol_type is a string -> char is needed to create title

    %Set all interpreters globally to Latex... what does groot mean?
    % set(groot,'TickLabelInterpreter','latex');
    % set(groot,'defaulttextinterpreter','latex');
    % set(groot,'defaultLegendInterpreter','latex');
    % set(gcf,'units','normalized');

    
    switch options.space

        case 'time'

            solget_options.eval = options.zaxis;                            % Tell solget what to evaluate
            solget_output = obj.solget(DYN,solget_options);                 % Get the solution(s)
            z              = solget_output.solution_eval;                   % Extract the fields from solget_output which are needed
            time           = solget_output.time;
            mu             = solget_output.mu;
            solget_options = solget_output.options;

            % Plot
            idx_mu = zeros(numel(solget_options.index)*size(z,2),2);        % Preallocate
            for k = 1:numel(solget_options.index)   % k loop: through the indices; j loop: needed, for the legend to work properly
                for j = 1:size(z,2)
                    counter = counter + 1;
                    idx_mu(counter,:) = [solget_options.index(k),round(mu(k)*100)/100];     % Round mu to 2 decimals
                    plot(time(:,1,k),z(:,j,k),'Color',obj.custom_color(options),'linestyle',linestyle);
                end
            end

            % Set title and labels (legend is set after switch...case because it is the same code for all cases, except quasi-periodic hypertime)
            title(append('Time Domain of ', upper(sol_type(1)), sol_type(2:end),' Solution'),'Interpreter','latex');
            xlabel('Time $t$','Interpreter','latex');
            xlim([min(time,[],'all'),max(time,[],'all')]);
            ylabel('State $z_i (t)$','Interpreter','latex');

            % Output
            output.z    = z;
            output.time = time;
            output.mu   = mu;


        case 'trajectory'                       % If case trajectroy is chosen: axis are always function_handles automatically

            solget_options.space = 'time';      % Trajectory is a time solution plotted in state space
            solget_options.eval = 'all';        % Tell solget what to evaluate (function handles apllied later)
            solget_output = obj.solget(DYN,solget_options);                     % Get the solution(s)
            s              = solget_output.solution_eval;                       % Extract the fields from solget_output which are needed
            mu             = solget_output.mu;
            solget_options = solget_output.options;
            idx_mu = zeros(numel(solget_options.index),2);                      % Preallocate
            x = zeros(size(s,1),1,numel(solget_options.index));                 % Preallocate
            z = zeros(size(s,1),1,numel(solget_options.index));                 % Preallocate

            % 3D Plot
            if isfield(options,'yaxis')                     

                % Plot
                y = zeros(size(s,1),1,numel(solget_options.index));             % Preallocate
                for k = 1:numel(solget_options.index)           % k loop: through the indices (no j loop needed since there is only one plot per solution)
                    x(:,1,k) = options.xaxis(s(:,:,k));         % The function_handle is applied at this stage - otherwise solget would need to be called 3 times
                    z(:,1,k) = options.zaxis(s(:,:,k));
                    y(:,1,k) = options.yaxis(s(:,:,k));
                    counter = counter + 1;
                    idx_mu(counter,:) = [solget_options.index(k),round(mu(k)*100)/100];   % Round mu to 2 decimals
                    plot3(x(:,1,k),z(:,1,k),y(:,1,k),'Color',obj.custom_color(options));
                end

                view(3);                % hold on also holds the view command... view(3) switches to 3D depiction
                zlabel('State $z_i$','Interpreter','latex');

                % Output
                output.z  = z;
                output.x  = x;
                output.y  = y;
                output.mu = mu;

            % 2D Plot
            else    
                
                % Plot
                for k = 1:numel(solget_options.index)       % k loop: through the indices (no j loop needed since there is only one plot per solution)
                    x(:,1,k) = options.xaxis(s(:,:,k));     % The function_handle is applied at this stage - otherwise solget would need to be called 2 times
                    z(:,1,k) = options.zaxis(s(:,:,k));
                    counter = counter + 1;
                    idx_mu(counter,:) = [solget_options.index(k),round(mu(k)*100)/100];     % Round mu to 2 decimals
                    plot(x(:,1,k),z(:,1,k),'Color',obj.custom_color(options),'linestyle',linestyle);
                end

                % Output
                output.z  = z;
                output.x  = x;                
                output.mu = mu;

            end

            % Set title and labels (legend is set after switch...case because it is the same code for all cases, except quasi-periodic hypertime)
            title(append('Trajectory of ', upper(sol_type(1)), sol_type(2:end),' Solution'),'Interpreter','latex');
            xlabel('State $z_j$','Interpreter','latex');
            ylabel('State $z_k$','Interpreter','latex');


        case 'hypertime'

            solget_options.eval = options.zaxis;                            % Tell solget what to evaluate
            solget_output = obj.solget(DYN,solget_options);                 % Get the solution(s)
            z              = solget_output.solution_eval;                   % Extract the fields from solget_output which are needed
            hypertime      = solget_output.hypertime;
            mu             = solget_output.mu;
            solget_options = solget_output.options;

            % Periodic solutions (solplot is not available for equilibrium solutions)
            if DYN.n_freq == 1
                
                % Plot
                idx_mu = zeros(numel(solget_options.index)*size(z,2),2);        % Preallocate
                for k = 1:numel(solget_options.index)       % k loop: through the indices; j loop: needed, for the legend to work properly
                    for j =1:size(z,2)
                        counter = counter + 1;
                        idx_mu(counter,:) = [solget_options.index(k),round(mu(k)*100)/100];   % Round mu to 2 decimals
                        plot(hypertime(:,1,k),z(:,j,k),'Color',obj.custom_color(options),'linestyle',linestyle);
                    end
                end

                % Set title and labels (legend is set after switch...case because it is the same code for all cases, except quasi-periodic hypertime)
                title(append('Hypertime Domain of ', upper(sol_type(1)), sol_type(2:end),' Solution'),'Interpreter','latex');
                xlabel('Hypertime $\theta$','Interpreter','latex')
                xlim([0,2*pi]);
                xticks([0,pi/2, pi,3/2*pi, 2.*pi])
                xticklabels({'$0$','$\frac{1}{2}\pi$','$\pi$','$\frac{3}{2}\pi$','$2 \pi$'})
                ylabel('State $z_i (\theta)$','Interpreter','latex');

            % Quasi-Periodic solutions
            elseif DYN.n_freq == 2

                for k = 1:numel(solget_options.index)       % k loop: through the indices; j loop: needed, for the legend to work properly
                    for j =1:size(z,3)

                        % Plot
                        counter = counter + 1;
                        if numel(solget_options.index)*size(z,3) > 1
                            subplot(numel(solget_options.index),size(z,3),counter);     % Create subplots if there are more than one solutions to be plotted
                        end
                        surf(hypertime(:,:,1,k),hypertime(:,:,2,k),z(:,:,j,k));
                        
                        % Set title, labels and ticks and view -> Set the information into the title for each subplot instead of a legend
                        set(gca,'TickLabelInterpreter','latex');    % Needs to be set again due to subplots
                        title({'Hypertime Domain',['IDX: ',num2str(solget_options.index(k)),...
                               '; $\mu_',num2str(DYN.act_param),'=',num2str(round(mu(k)*100)/100),'$']},'Interpreter','latex');
                        xlim([0,2*pi]);
                        xticks([0,pi/2, pi,3/2*pi, 2.*pi]);
                        xticklabels({'$0$','$\frac{1}{2}\pi$','$\pi$','$\frac{3}{2}\pi$','$2\pi$'});
                        ylim([0,2*pi]);
                        yticks([0,pi/2, pi,3/2*pi, 2.*pi]);
                        yticklabels({'$0$','$\frac{1}{2}\pi$','$\pi$','$\frac{3}{2}\pi$','$2 \pi$'});
                        if numel(solget_options.index)*size(z,3) == 1
                            xlabel('Hypertime $\theta_1$','Interpreter','latex')
                            ylabel('Hypertime $\theta_2$','Interpreter','latex')
                        else
                            xlabel('$\theta_1$','Interpreter','latex')          % Omit 'Hypertime' in axes label if there are more than 1 plots
                            ylabel('$\theta_2$','Interpreter','latex')          % Omit 'Hypertime' in axes label if there are more than 1 plots
                        end
                        zlabel('$Z_i (\mbox{\boldmath$\theta$})$','Interpreter','latex');
                        view([-30,50]);     

                    end
                end

            end

            % Output
            output.z         = z;
            output.hypertime = hypertime;
            output.mu        = mu;


        case 'frequency'

            solget_options.eval = options.zaxis;                            % Tell solget what to evaluate
            solget_output = obj.solget(DYN,solget_options);                 % Get the solution(s)
            z              = solget_output.amplitude;                       % Extract the fields from solget_output which are needed
            alpha          = solget_output.angle;
            f              = solget_output.frequency;
            mu             = solget_output.mu;
            solget_options = solget_output.options;

            % Plot
            idx_mu = zeros(numel(solget_options.index)*size(z,2),2);        % Preallocate
            for k = 1:numel(solget_options.index)       % k loop: through the indices; j loop: needed, for the legend to work properly
                for j = 1:size(z,2)
                    counter = counter + 1;
                    idx_mu(counter,:) = [solget_options.index(k),round(mu(k)*100)/100];   % Round mu to 2 decimal places
                    plot(f(:,1,k),z(:,j,k),'Color',obj.custom_color(options),'linestyle',linestyle);
                end
            end

            % Set title and labels (legend is set after switch...case because it is the same code for all cases, except quasi-periodic hypertime)
            title(append('Frequency Domain of ', upper(sol_type(1)), sol_type(2:end),' Solution'),'Interpreter','latex');
            xlabel('Angular Frequency $\omega$','Interpreter','latex');
            ylabel('Absolute Amplitude $|\mathcal{F}_i (\omega)|$','Interpreter','latex');

            % Output -> TODO: Implement for 3D
            output.amplitude = z;
            output.angle     = alpha;
            output.frequency = f;
            output.mu        = mu;

    end


    % Set the legend (except for quasi-periodic solutions in hypertime space)
    if ~(strcmpi(options.space,'hypertime') && DYN.n_freq == 2)
        LegStr = cell(1,size(idx_mu,1));                            % Allocate cell memory
        for k = 1:size(idx_mu,1)
            LegStr{1,k} = ['IDX: ',num2str(idx_mu(k,1)),'; $\mu_',num2str(DYN.act_param),'=',num2str(idx_mu(k,2)),'$'];
        end
        legend([LegStr_old,LegStr],'Interpreter','latex');          % Set the legend
    end


end