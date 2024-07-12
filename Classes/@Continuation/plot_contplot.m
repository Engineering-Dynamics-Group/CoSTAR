    %This is a method of the continuation class and plots a continuation plot
    %during a continuation
    %
    %@obj:  Continuation class object
    %@S:    Soltuion subclass object
    %@DYN:  DynamicalSystem class object
    
function plot_contplot(obj,S,DYN)
    
    % Adapt the resolution depending on the solution type (resolution is listed w.r.t. torus coordinate)
    if strcmpi(DYN.sol_type,'quasiperiodic')
        resolution = 50;
    else
        resolution = 200;
    end
    

    % Select color
    mycolor = 'b';                                                          % Stability was not computed or stable solution
    if (strcmpi(DYN.stability,'on') == 1) && (obj.p_n_unstable_1 > 0)       % Unstable solution
        mycolor = 'r';
    end
    % if (strcmpi(DYN.stability,'on') == 1) && (obj.p_stability_flag == 2)    % Stability cannot be determined based on multipliers
    %     mycolor = S.plot_color.magenta;
    % end
    max_idx = numel(S.mu);


    % Start plotting after first continuation loop
    if obj.p_local_cont_counter == 1
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'color',mycolor,'resolution',resolution);
        [s,mu] = S.contplot(DYN,opts);
        obj.p_axes_values_old = [mu(2),s(2)];       % Save axes values of second solution (i.e. first solution after initial solution) in order to use them in the next loop
        tmp = ylim;
        obj.p_limit = 2*tmp(2);


    % Plotting for: No change in stability OR stability is not computed
    elseif (obj.p_n_unstable_1 == obj.p_n_unstable_0) || strcmpi(DYN.stability,'off')  
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'figure',gcf,'color',mycolor,'resolution',resolution,'axes_values_old',obj.p_axes_values_old);
        [s,mu] = S.contplot(DYN,opts);
        obj.p_axes_values_old = [mu(2),s(2)];       % Save axes values of second solution in order to use them in the next loop


    % Plotting for: Change in stability, but the bifurcation points were NOT iterated
    % Plot in dark grey if stability changes from stable to unstable or vice versa (The stability behaviour is unclear in this section ...
    % since the location of the BFP is unknown -> plotting in blue or red would be misleading)
    elseif isfield(DYN.opt_stability,'iterate_bfp') && strcmpi(DYN.opt_stability.iterate_bfp,'off') && (obj.p_n_unstable_1 ~= obj.p_n_unstable_0)
        if (obj.p_n_unstable_1 == 0) || (obj.p_n_unstable_0 == 0)
            mycolor = S.plot_color.dark_grey;
        end
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'figure',gcf,'color',mycolor,'resolution',resolution,'axes_values_old',obj.p_axes_values_old);
        [s,mu] = S.contplot(DYN,opts);
        obj.p_axes_values_old = [mu(2),s(2)];       % Save axes values of second solution in order to use them in the next loop
    

    % Plotting for: Change in stability AND the bifurcation points were iterated
    % The following three cases are necessary, since 2 line segments (3 indices) need to be plotted, 
    % if a bifurcation point occurred in between two curve points. The different cases are needed for the correct color.
    elseif (obj.p_n_unstable_1 > obj.p_n_unstable_0) && obj.p_n_unstable_0==0
        opts = costaropts('zaxis','max2','index',double([max_idx-2,max_idx-1]),'figure',gcf,'color','b','resolution',resolution,'axes_values_old',obj.p_axes_values_old);
        [s,mu] = S.contplot(DYN,opts);
        obj.p_axes_values_old = [mu(2),s(2)];       % Save axes values of second solution in order to use them in the next contplot call
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'figure',gcf,'color','r','resolution',resolution,'axes_values_old',obj.p_axes_values_old);
        [s,mu] = S.contplot(DYN,opts);
        obj.p_axes_values_old = [mu(2),s(2)];       % Save axes values of second solution in order to use them in the next loop
        
    elseif (obj.p_n_unstable_1 < obj.p_n_unstable_0) && obj.p_n_unstable_1==0
        opts = costaropts('zaxis','max2','index',double([max_idx-2,max_idx-1]),'figure',gcf,'color','r','resolution',resolution,'axes_values_old',obj.p_axes_values_old);
        [s,mu] = S.contplot(DYN,opts);
        obj.p_axes_values_old = [mu(2),s(2)];       % Save axes values of second solution in order to use them in the next contplot call
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'figure',gcf,'color','b','resolution',resolution,'axes_values_old',obj.p_axes_values_old);
        [s,mu] = S.contplot(DYN,opts);
        obj.p_axes_values_old = [mu(2),s(2)];       % Save axes values of second solution in order to use them in the next loop
    
    elseif (obj.p_n_unstable_1 ~= obj.p_n_unstable_0)
        opts = costaropts('zaxis','max2','index',double([max_idx-2,max_idx-1]),'figure',gcf,'color',mycolor,'resolution',resolution,'axes_values_old',obj.p_axes_values_old);
        [s,mu] = S.contplot(DYN,opts);
        obj.p_axes_values_old = [mu(2),s(2)];       % Save axes values of second solution in order to use them in the next contplot call
        opts = costaropts('zaxis','max2','index',double([max_idx-1,max_idx]),'figure',gcf,'color',mycolor,'resolution',resolution,'axes_values_old',obj.p_axes_values_old);
        [s,mu] = S.contplot(DYN,opts);
        obj.p_axes_values_old = [mu(2),s(2)];       % Save axes values of second solution in order to use them in the next loop

    end
    

    % Set y limits (x limits have already been set in contplot)
    tmp = [obj.p_limit,0.9.*min(s),1.1.*max(s)];
    obj.p_limit = max(tmp);
    ylim([0,obj.p_limit])

       
end