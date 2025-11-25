% This is a method of the continuation class and plots a continuation plot during a continuation
%
% @obj:  Continuation class object
% @S:    Solution subclass object
% @DYN:  DynamicalSystem class object

function plot_contplot(obj,S,DYN)

    % Select color
    if ~isempty(obj.p_stability_flag) && ((obj.p_stability_flag == 0) || (obj.p_stability_flag_old == 0))   
        mycolor = S.plot_color.orange;                                  % If stability computation failed for latest or previous solution: plot in orange
    elseif strcmpi(DYN.stability,'on') && (obj.p_n_unstable_1 == 0)     % Stable solution
        mycolor = 'g';
    elseif strcmpi(DYN.stability,'on') && (obj.p_n_unstable_1 > 0)      % Unstable solution
        mycolor = 'r';
    else
        mycolor = 'b';                                                  % Stability was not computed
    end
    % if (strcmpi(DYN.stability,'on') == 1) && (obj.p_stability_flag == 2)    % Stability cannot be determined based on multipliers
    %     mycolor = S.plot_color.magenta;
    % end
    max_idx = numel(S.mu);


    % Start plotting after first continuation loop
    if obj.p_local_cont_counter == 1
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'color',mycolor);
        contplot_output = S.contplot(DYN,opts);
        obj.p_axes_values_old = [contplot_output.mu(2),contplot_output.z(2)];   % Save axes values of second solution in order to use them in the next loop
        tmp = ylim;
        obj.p_limit = 2*tmp(2);


    % Plotting for: No change in stability OR stability is not computed
    elseif (obj.p_n_unstable_1 == obj.p_n_unstable_0) || strcmpi(DYN.stability,'off')
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'figure',gcf,'color',mycolor,'axes_values_old',obj.p_axes_values_old);
        contplot_output = S.contplot(DYN,opts);
        obj.p_axes_values_old = [contplot_output.mu(2),contplot_output.z(2)];   % Save axes values of second solution in order to use them in the next loop


    % Plotting for: Stability computation failed for latest or previous solution
    elseif (obj.p_stability_flag == 0) || (obj.p_stability_flag_old == 0)
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'figure',gcf,'color',S.plot_color.orange,'axes_values_old',obj.p_axes_values_old);
        contplot_output = S.contplot(DYN,opts);
        obj.p_axes_values_old = [contplot_output.mu(2),contplot_output.z(2)];   % Save axes values of second solution in order to use them in the next loop


    % Plotting for: Change in stability, but the bifurcation points were NOT iterated OR iteration failed
    % Plot in pale grey if stability changes from stable to unstable or vice versa (The stability behaviour is unclear in this section ...
    % since the location of the BFP is unknown -> plotting in blue or red would be misleading)
    elseif (obj.p_n_unstable_1 ~= obj.p_n_unstable_0) && (isempty(obj.p_newton_flag_bfp) || obj.p_newton_flag_bfp < 1)
        if (obj.p_n_unstable_1 == 0) || (obj.p_n_unstable_0 == 0)
            mycolor = S.plot_color.pale_grey;
        end
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'figure',gcf,'color',mycolor,'axes_values_old',obj.p_axes_values_old);
        contplot_output = S.contplot(DYN,opts);
        obj.p_axes_values_old = [contplot_output.mu(2),contplot_output.z(2)];   % Save axes values of second solution in order to use them in the next loop
    

    % Plotting for: Change in stability AND the bifurcation points were successfully iterated
    % The following three cases are necessary, since 2 line segments (3 indices) need to be plotted, 
    % if a bifurcation point occurred in between two curve points. The different cases are needed for the correct color.
    elseif (obj.p_n_unstable_1 > obj.p_n_unstable_0) && obj.p_n_unstable_0==0
        opts = costaropts('zaxis','max2','index',double([max_idx-2,max_idx-1]),'figure',gcf,'color','g','axes_values_old',obj.p_axes_values_old);
        contplot_output = S.contplot(DYN,opts);
        obj.p_axes_values_old = [contplot_output.mu(2),contplot_output.z(2)];   % Save axes values of second solution in order to use them in the next contplot call
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'figure',gcf,'color','r','axes_values_old',obj.p_axes_values_old);
        contplot_output = S.contplot(DYN,opts);
        obj.p_axes_values_old = [contplot_output.mu(2),contplot_output.z(2)];   % Save axes values of second solution in order to use them in the next loop
        
    elseif (obj.p_n_unstable_1 < obj.p_n_unstable_0) && obj.p_n_unstable_1==0
        opts = costaropts('zaxis','max2','index',double([max_idx-2,max_idx-1]),'figure',gcf,'color','r','axes_values_old',obj.p_axes_values_old);
        contplot_output = S.contplot(DYN,opts);
        obj.p_axes_values_old = [contplot_output.mu(2),contplot_output.z(2)];   % Save axes values of second solution in order to use them in the next contplot call
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'figure',gcf,'color','g','axes_values_old',obj.p_axes_values_old);
        contplot_output = S.contplot(DYN,opts);
        obj.p_axes_values_old = [contplot_output.mu(2),contplot_output.z(2)];   % Save axes values of second solution in order to use them in the next loop
    
    elseif (obj.p_n_unstable_1 ~= obj.p_n_unstable_0)
        opts = costaropts('zaxis','max2','index',double([max_idx-2,max_idx-1]),'figure',gcf,'color',mycolor,'axes_values_old',obj.p_axes_values_old);
        contplot_output = S.contplot(DYN,opts);
        obj.p_axes_values_old = [contplot_output.mu(2),contplot_output.z(2)];   % Save axes values of second solution in order to use them in the next contplot call
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'figure',gcf,'color',mycolor,'axes_values_old',obj.p_axes_values_old);
        contplot_output = S.contplot(DYN,opts);
        obj.p_axes_values_old = [contplot_output.mu(2),contplot_output.z(2)];   % Save axes values of second solution in order to use them in the next loop

    end
    

    % Set y limits (x limits have already been set in contplot)
    tmp = [obj.p_limit,0.9.*min(contplot_output.z),1.1.*max(contplot_output.z)];
    obj.p_limit = max(tmp);
    ylim([0,obj.p_limit])

    drawnow
       
    drawnow
    
end