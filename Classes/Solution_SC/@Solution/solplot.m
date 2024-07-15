% This function plots prescribed solution points in a 2D or 3D plot
% necessary input is:
% @obj                 Object of Solution-Class
% @DYN                 Object of DynamicalSystem-Class
% @options             struct containing options for plotting/fieldnames:
%
%   [s_out,varargout] = solplot(obj,DYN,options)
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
%                       'index'         1xn array indiciating the indices
%                                       of the solution curve, where the
%                                       solution is evaluated.
%                                       Default: All solutions are given
%                                       back
%                       'mu'            1xn array indiciaing the mu values
%                                       of the solution curve, where the
%                                       solution is evaluated. If the exact
%                                       value is not present, the closest
%                                       value is used. If these leads to
%                                       doubled values, they are ignored.
%                                       Mu might lead to ambuguities in
%                                       case of overhanging curves
%                       'color'         Define a color for the plot either by 
%                                       providing an rgb array or a char 
%                                       'r','g','b','c','m','y' or 'k'
% Output Arguments:
% @s_out:               Solution vector in the dimension [options.resolution,
%                       statespace dimension, number of curve points].
%                       statespace dimension depends on options.eval
% @varargout:           Contains the data, which is plotted
%                       options.space = 'time'
%                           varargout{1,1}: time
%                           varargout{1,2}: zaxis
%                           varargout{1,3}: mu
%                           varargout{1,4}: empty
%                        options.space = 'trajectory' 2D
%                           varargout{1,1}: xaxis
%                           varargout{1,2}: zaxis
%                           varargout{1,3}: mu
%                           varargout{1,4}: empty
%                        options.space = 'trajectory' 3D
%                           varargout{1,1}: xaxis
%                           varargout{1,2}: zaxis
%                           varargout{1,3}: yaxis
%                           varargout{1,4}: mu
%                        options.space = 'hypertime' 2D
%                           varargout{1,1}: hypertime
%                           varargout{1,2}: zaxis
%                           varargout{1,3}: mu
%                           varargout{1,4}: empty
%                        options.space = 'hypertime' 3D (not implemented yet)
%                           varargout{1,1}: hypertime 1
%                           varargout{1,2}: hypertime 1
%                           varargout{1,3}: zaxis
%                           varargout{1,4}: mu
%                        options.space = 'frequency' 2D
%                           varargout{1,1}: frequency
%                           varargout{1,2}: zaxis
%                           varargout{1,3}: mu
%                           varargout{1,4}: empty
%                        options.space = 'frequency' 3D (not implemented yet)
%                           varargout{1,1}: frequency along basefrequency 1
%                           varargout{1,2}: frequency along basefrequency 2
%                           varargout{1,3}: zaxis
%                           varargout{1,4}: mu

%Example:   options = struct('xaxis',@(z)z(:,1),'zaxis',@(z)z(:,2),'space','trajectory','mu',[0.5,1,1.5]);
%           [x,y] = S.soplot(DYN,options);


function varargout = solplot(obj,DYN,options)

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

    hold on;                                        % Allow multiple lines to be plotted and enable grid
    counter = 0;                                    % Counter needed for legend
    %Set all interpreters globally to Latex... what does groot mean?
    % set(groot,'TickLabelInterpreter','latex');
    % set(groot,'defaulttextinterpreter','latex');
    % set(groot,'defaultLegendInterpreter','latex');
    % set(gcf,'units','normalized');

    
    switch options.space

        case 'time'

            solget_options.eval = options.zaxis;                            % Tell solget what to evaluate
            [z,mu,time,solget_options] = obj.solget(DYN,solget_options);    % Get the solution(s)

            % Plot
            LegStr_old = cell(get(findobj(gcf,'type','Legend'),'String'));  % Get existing legend
            idx_mu = zeros(numel(solget_options.index)*size(z,2),2);        % Preallocate
            for k = 1:numel(solget_options.index)   % k loop: through the indices; j loop: needed, for the legend to work properly
                for j = 1:size(z,2)
                    counter = counter + 1;
                    idx_mu(counter,:) = [solget_options.index(k),round(mu(k)*100)/100];     % Round mu to 2 decimals
                    plot(time(:,1,k),z(:,j,k),'Color',obj.custom_color(options),'linestyle',linestyle);
                end
            end

            % Set title, labels, limits and legend
            title('Time Domain','Interpreter','latex');
            xlabel('Time $t$','Interpreter','latex');
            xlim([min(time,[],'all'),max(time,[],'all')]);
            ylabel('State $z_i (t)$','Interpreter','latex');
            LegStr = cell(1,size(idx_mu,1));                            % Allocate cell memory
            for k = 1:size(idx_mu,1)
                LegStr{1,k} = ['IDX: ',num2str(idx_mu(k,1)),'; $\mu_',num2str(DYN.act_param),'=',num2str(idx_mu(k,2)),'$'];
            end
            legend([LegStr_old,LegStr],'Interpreter','latex');          % Set the legend

            % Output
            varargout = cell(1,4);
            varargout{1,1} = time;
            varargout{1,2} = z;
            varargout{1,3} = mu;
            varargout{1,4} = [];


        case 'trajectory'                       % If case trajectroy is chosen: axis are always function_handles automatically

            solget_options.space = 'time';      % Trajectory is a time solution plotted in state space
            solget_options.eval = 'all';        % Tell solget what to evaluate (function handles apllied later)
            [s_traj,mu,~,solget_options] = obj.solget(DYN,solget_options);      % Get the solution(s)
            LegStr_old = cell(get(findobj(gcf,'type','Legend'),'String'));      % Get existing legend

            % 3D Plot
            if isfield(options,'yaxis')                     

                % Plot
                idx_mu = zeros(numel(solget_options.index),2);            % Preallocate
                for k = 1:numel(solget_options.index)       % k loop: through the indices (no j loop needed since there is only one plot per solution)
                    x = options.xaxis(s_traj(:,:,k));       % The function_handle is applied at this stage - otherwise solget would need to be called 3 times
                    z = options.zaxis(s_traj(:,:,k));
                    y = options.yaxis(s_traj(:,:,k));
                    counter = counter + 1;
                    idx_mu(counter,:) = [solget_options.index(k),round(mu(k)*100)/100];   % Round mu to 2 decimals
                    plot3(x,z,y,'Color',obj.custom_color(options));
                end

                view(3);                % hold on also holds the view command... view(3) switches to 3D depiction
                zlabel('State $z_i$','Interpreter','latex');

                % Output
                varargout = cell(1,4);
                varargout{1,1} = x;
                varargout{1,2} = z;
                varargout{1,3} = y;
                varargout{1,4} = mu;

            % 2D Plot
            else    
                
                % Plot
                idx_mu = zeros(numel(solget_options.index),2);                      % Preallocate
                for k = 1:numel(solget_options.index)       % k loop: through the indices (no j loop needed since there is only one plot per solution)
                    x = options.xaxis(s_traj(:,:,k));       % The function_handle is applied at this stage - otherwise solget would need to be called 2 times
                    z = options.zaxis(s_traj(:,:,k));
                    counter = counter + 1;
                    idx_mu(counter,:) = [solget_options.index(k),round(mu(k)*100)/100];     % Round mu to 2 decimals
                    plot(x,z,'Color',obj.custom_color(options),'linestyle',linestyle);
                end

                % Output
                varargout = cell(1,4);
                varargout{1,1} = x;
                varargout{1,2} = z;
                varargout{1,3} = mu;
                varargout{1,4} = [];

            end

            % Set title, labels and legend
            title('Trajectory','Interpreter','latex');
            xlabel('State $z_j$','Interpreter','latex');
            ylabel('State $z_k$','Interpreter','latex');
            LegStr = cell(1,size(idx_mu,1));                            % Allocate cell memory
            for k = 1:size(idx_mu,1)
                LegStr{1,k} = ['IDX: ',num2str(idx_mu(k,1)),'; $\mu_',num2str(DYN.act_param),'=',num2str(idx_mu(k,2)),'$'];
            end
            legend([LegStr_old,LegStr],'Interpreter','latex');          % Set the legend


        case 'hypertime'

            solget_options.eval = options.zaxis;                                % Tell solget what to evaluate
            [z,mu,hypertime,solget_options] = obj.solget(DYN,solget_options);   % Get the solution(s)

            % Periodic solutions (solplot is not available for equilibrium solutions)
            if DYN.n_freq == 1

                % Plot
                LegStr_old = cell(get(findobj(gcf,'type','Legend'),'String'));  % Get existing legend
                idx_mu = zeros(numel(solget_options.index)*size(z,2),2);        % Preallocate
                for k = 1:numel(solget_options.index)       % k loop: through the indices; j loop: needed, for the legend to work properly
                    for j =1:size(z,2)
                        counter = counter + 1;
                        idx_mu(counter,:) = [solget_options.index(k),round(mu(k)*100)/100];   % Round mu to 2 decimals
                        plot(hypertime(:,1,k),z(:,j,k),'Color',obj.custom_color(options),'linestyle',linestyle);
                    end
                end

                % Set title, labels, ticks and legend
                title('Hypertime Domain','Interpreter','latex');
                xlabel('Hypertime $\theta$','Interpreter','latex')
                xlim([0,2*pi]);
                xticks([0,pi/2, pi,3/2*pi, 2.*pi])
                xticklabels({'$0$','$\frac{1}{2}\pi$','$\pi$','$\frac{3}{2}\pi$','$2 \pi$'})
                ylabel('State $z_i (\theta)$','Interpreter','latex');
                LegStr = cell(1,size(idx_mu,1));                            % Allocate cell memory
                for k = 1:size(idx_mu,1)
                    LegStr{1,k} = ['IDX: ',num2str(idx_mu(k,1)),'; $\mu_',num2str(DYN.act_param),'=',num2str(idx_mu(k,2)),'$'];
                end
                legend([LegStr_old,LegStr],'Interpreter','latex');          % Set the legend

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
                        zlabel('$z_i (\theta)$','Interpreter','latex');
                        view([-30,50]);     

                    end
                end

            end

            % Output
            varargout = cell(1,4);
            varargout{1,1} = hypertime;
            varargout{1,2} = z;
            varargout{1,3} = mu;
            varargout{1,4} = [];


        case 'frequency'

            solget_options.eval = options.zaxis;                            % Tell solget what to evaluate
            [z,mu,f,solget_options] = obj.solget(DYN,solget_options);       % Get the solution(s)

            % Plot
            LegStr_old = cell(get(findobj(gcf,'type','Legend'),'String'));  % Get existing legend
            idx_mu = zeros(numel(solget_options.index)*size(z,2),2);        % Preallocate
            for k = 1:numel(solget_options.index)       % k loop: through the indices; j loop: needed, for the legend to work properly
                for j = 1:size(z,2)
                    counter = counter + 1;
                    idx_mu(counter,:) = [solget_options.index(k),round(mu(k)*100)/100];   % Round mu to 2 decimal places
                    plot(f(:,1,k),z(:,j,k),'Color',obj.custom_color(options),'linestyle',linestyle);
                end
            end

            % Set title, labels and legend
            title('Frequency Domain','Interpreter','latex');
            xlabel('Angular Frequency $\omega$','Interpreter','latex');
            ylabel('Absolute Amplitude $|\mathcal{F}(z_i)|$','Interpreter','latex');
            LegStr = cell(1,size(idx_mu,1));                            % Allocate cell memory
            for k = 1:size(idx_mu,1)
                LegStr{1,k} = ['IDX: ',num2str(idx_mu(k,1)),'; $\mu_',num2str(DYN.act_param),'=',num2str(idx_mu(k,2)),'$'];
            end
            legend([LegStr_old,LegStr],'Interpreter','latex');          % Set the legend

            % Output -> TODO: Implement for 3D
            varargout = cell(1,4);
            varargout{1,1} = f;
            varargout{1,2} = z;
            varargout{1,3} = mu;
            varargout{1,4} = [];

    end


end