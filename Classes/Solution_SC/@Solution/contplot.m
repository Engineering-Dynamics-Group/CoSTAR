% This function plots a  solution curve in a 2D plot. It evaluates all
% solution curve points contained in SOLUTION subclass object based on the
% hypertime depiction of the solution object
%
% Input Arguments
% @obj:                 Object of Solution class
% @DYN:                 Object of DynamicalSystem class
% @options:             struct contain options for calculating solution
%
%                       Mandatory fields:
%                       'zaxis':    -->'function_handle:' User supplied
%                                       function handle. Example: options('zaxis',@(z)max(z(:,1))
%                                       function handle must give back a
%                                       scalar.
%                                   --> 'min2': plots the
%                                       minimum of the euclidean norm,
%                                       which is evaluated along the states
%                                       of an output vector over the
%                                       continuation parameter mu
%                                   --> 'mean2': plots the
%                                       mean of the euclidean norm,
%                                       which is evaluated along the states
%                                       of an output vector over the
%                                       continuation parameter mu
%                                   --> 'max2': gives back the
%                                       maximum of the euclidean norm,
%                                       which is evaluated along the states
%                                       of an output vector over the
%                                       continuation parameter mu
%
%                       Optional fields:
%                       'resolution'    1x1 or 1xnumb_base_frq double array
%                                       defining the number of points in
%                                       time or hypertime at which the
%                                       solution is evaluated. This is the
%                                       total number of points (see 'interval')
%                                       Default value: 200
%                                       This might increase the precision
%                                       of the continuation curve
%
%                       'figure'        handle to a existing figure. If
%                                       supplied, solplot plots into the
%                                       existing figure environment.
%                       'color'         Define a color for the plot either by 
%                                       providing an rgb array or a char 
%                                       'r','g','b','c','m','y' or 'k'
%
% Output Arguments:
%
% @varargout:           Contains the data, which is plotted
%                           varargout{1,1}: zaxis
%                           varargout{1,2}: mu

%Example:   options = struct('zaxis',@(z)z(:,1));
%           [x,y] = S.contplot(DYN,options);


function varargout = contplot(obj,DYN,options)

%% Gatekeeper check of input: Only the input not checked by solget_gatekeeper is checked here
if ~isfield(options,'axes_values_old')                  % When this field is present, contplot is called from the plot_contplot method during a continuation
    options = obj.contplot_gatekeeper(DYN,options);     % The gatekeeper can be skipped in this case since the developers make sure that the options struct is fine
end

varargout = cell(1,2);                          % Needed

% Add the complete index array, if none or 'all' is specified (this is easier for the plotting and the bifurcations below)
if ~isfield(options,'index') || strcmpi(options.index,'all')
    options.index = 1:numel(obj.mu); 
end           



%% Set the solget options structure
solget_options = options;                       % solget_options is a subset of solplot options.
solget_options.space = 'hypertime';


% Remove fields that are not allowed in solget options structure
solget_options = rmfield(solget_options,'zaxis');   

if isfield(solget_options,'figure'); solget_options = rmfield(solget_options,'figure'); end
if isfield(solget_options,'color'); solget_options = rmfield(solget_options,'color'); end
if isfield(solget_options,'linestyle'); solget_options = rmfield(solget_options,'linestyle'); end


% This statement is going to work, since all inputs have been checked in the gatekeeper
if isa(options.zaxis,'char')                    % This is the case for min2, mean2, max2
    solget_options.eval = 'euclidean';
else
    solget_options.eval = options.zaxis;
end



%% Get the data
if isfield(options,'axes_values_old')                   % When this field is present, contplot is called from the plot_contplot method during a continuation
    solget_options.index = options.index(2);            % First index can be skipped since this point has already been evaluated 
    [s,mu] = obj.solget(DYN,solget_options);            % Only the new solution has to be evaluated
else
    [s,mu] = obj.solget(DYN,solget_options);            % contplot is called by the user
end


if isa(options.zaxis,'function_handle')         % In this case, contplot_gatekeeper has checked the correct output dimension of the function handle
    
    if DYN.n_freq == 0
        s_out = s; 
    elseif DYN.n_freq == 1
        s_out = permute(s,[2,3,1]);             
    else
        s_out = permute(s,[3,4,1,2]);
    end

    TtlStr = 'User Function Handle Applied to Hypertime Manifold';
    ystr = 'f(z)';


else
    switch options.zaxis     

        case 'min2'

            if DYN.n_freq  == 0 
                s_out = s;                                  % There is no minimum value for a point in state space
            elseif DYN.n_freq ==1 
                s_out = permute(min(s,[],1),[2,3,1]);       % Minimum on hypertime domain [0,2*pi]. Third array dim of s indicates solution point (index)
            else
                s_out = min(reshape(s,size(s,1)*size(s,2),size(s,4)),[],1);     % Minimum on hypertime domain [0,2*pi]^2. Fourth array dim of s indicates solution point (index)
            end

            TtlStr = 'Minimum of Euclidean Norm of Hypertime Manifold';
            ystr = 'min(norm(z,2))';
            

        case 'mean2'

            if DYN.n_freq  == 0 
                s_out = s;                                  % There is no mean value for a point in state space
            elseif DYN.n_freq == 1
                s_out = permute(mean(s,1),[2,3,1]);         % Mean on hypertime domain [0,2*pi]. Third array dim of s indicates solution point (index)
            else
                s_out = mean(reshape(s,size(s,1)*size(s,2),size(s,4)),1);    % Mean on hypertime domain [0,2*pi]^2. Fourth array dim of s indicates solution point (index)
            end

            TtlStr = 'Mean Value of Euclidean Norm of Hypertime Manifold';
            ystr = 'mean(norm(z,2))';

            
        case 'max2'
                
            if DYN.n_freq  == 0
                s_out = s;                                  % There is no maximum value for a point in state space
            elseif DYN.n_freq == 1
                s_out = permute(max(s,[],1),[2,3,1]);       % Maximum on hypertime domain [0,2*pi]. Third array dim of s indicates solution point (index)
            else
                s_out = max(reshape(s,size(s,1)*size(s,2),size(s,4)),[],1);     % Maximum on hypertime domain [0,2*pi]^2. Fourth array dim of s indicates solution point (index)
            end

            TtlStr = 'Maximum of Euclidean Norm of Hypertime Manifold';
            ystr = 'max(norm(z,2))';
            
    end
    
end


if isfield(options,'axes_values_old')               % When this field is present: mu and s_out are scalars (latest solution), but we need them to be a [1x2] vector of the latest two solutions for the plot
    mu = [options.axes_values_old(1),mu];           % Build the correct mu vector for the plot
    s_out = [options.axes_values_old(2),s_out];     % Build the correct s_out vector for the plot
end



%% Look for bifurcation within index range and get the bifurcation indices contained within

if ~isempty(obj.bifurcation)

    bfp_bool = ismember(obj.bifurcation.index,options.index);       % See if the bifurcation indices are within the requested indices
    bfp_idx = obj.bifurcation.index(logical(bfp_bool));             % Get the bifurcation points which are within the requested indices

    % Get the local index of the bifurcation point (bfp_idx describes the index with respect to the complete continuation)
    bifurc_idx_local = zeros(1,numel(bfp_idx));                     % Initialize
    for k = 1:numel(bfp_idx)
        bifurc_idx_local(k) = find(options.index==bfp_idx(k));      % bifurc_idx_local describes the bfp index with respect to [mu,s_out], i.e. the requested indeces to plot
    end

else

    bifurc_idx_local = [];      % Create an empty array (needed below)

end



%% Plotting

% Settings
if isfield(options,'figure')
    figure(options.figure);                     % Use existing figure
    set(gcf,'DefaultLineLineWidth',2)
else
    figure;                                     % Create new plot figure
    set(gcf,'Color','w','NumberTitle','off','DefaultLineLineWidth',2,'Name','Plot by SOLUTION.contplot');
    set(gca,'TickLabelInterpreter','latex');    % Set tick label interpreter to LaTeX
    grid on;
end

hold on;                                        % Allow multiple lines to be plotted and enable grid

if isfield(options,'linestyle')                 % Define the linestyle
    linestyle = options.linestyle; 
else
    linestyle = '-';
end


% Plot
if isfield(options,'color') || numel(obj.n_unstable)==0                             % Plot if custom color was set OR stability was not calculated

    plot(mu,s_out,'color',obj.custom_color(options),'linestyle',linestyle); 


elseif (numel(obj.n_unstable)~=0) && isfield(DYN.opt_stability,'iterate_bfp') && strcmpi(DYN.opt_stability.iterate_bfp,'off')   % Plot if stability was calculated, but bifurcation points were NOT iterated

    % Create the idx vector, which defines the start and end indices of the different parts (regarding the stability behaviour) of the solution path
    % idx_n_unstable_change contains the indices before and after the number of unstable multipliers change (needed for idx)
    n_unstable = obj.n_unstable(1,options.index);       % Get the number of unstable multipliers
    idx_diff_n_unstable = find(diff(n_unstable));       % Find the locations where number of unstable multipliers change
    idx_n_unstable_change = sort([idx_diff_n_unstable, idx_diff_n_unstable+ones(1,numel(idx_diff_n_unstable))],'ascend');
    idx = unique([1 idx_n_unstable_change numel(options.index)]);       % unique() is needed in case that idx_n_unstable_change contains the values 1 or numel(options.index) 

    % Loop for the parts of the solution path (regarding the stability behaviour)
    for k = 1:(numel(idx)-1)

        % Set the color depending on the stability behaviour
        if n_unstable(idx(k)) == 0;     options.color = 'b';     else;     options.color = 'r';     end

        % Security check: Are two idx points following each other? This is the case if stability behaviour changes
        if (idx(k+1) - idx(k)) == 1
            if (n_unstable(idx(k)) == 0) || (n_unstable(idx(k+1)) == 0)     % If stability changes from unstable to stable or vice versa, ...
                options.color = obj.plot_color.dark_grey;                   % ... a bifurcation occured but the exact bifurcation point is unknown (not iterated). Therefore: Plot in dark grey
            end
        end

        plot(mu(idx(k):idx(k+1)),s_out(idx(k):idx(k+1)),'color',obj.custom_color(options),'linestyle',linestyle);

    end

else        % Plot if stability was calculated AND the bifurcation points were iterated

    % Create the idx vector, which defines the start and end indices of the different parts (regarding the stability behaviour) of the solution path
    n_unstable = obj.n_unstable(1,options.index);               % Get the number of unstable multipliers
    idx = unique([1, bifurc_idx_local, numel(options.index)]);  % unique() is needed in case that first or last requested index corresponds to a bfp
    if isempty(bifurc_idx_local)                                % If bifurc_idx_local is empty ...
        bifurc_idx_local = 0;                                   % ... set it to 0 because code below does not work if bifurc_idx_local is empty
    end

    % Loop for the parts of the solution path (regarding the stability behaviour)
    for k = 1:(numel(idx)-1)                    

        % Security check: Are two idx points following each other?
        if ((idx(k+1) - idx(k)) == 1)            
            if (k == 1) && (bifurc_idx_local(1) ~= 1)           % First section to be plotted AND section does not start at a BFP
                idx_stab = 1;                                   % Stability behaviour can be determined at the first curve point
            elseif (k == numel(idx)-1) && (bifurc_idx_local(end) ~= idx(end))   % Last section to be plotted AND section does not end at a BFP
                idx_stab = idx(end);                            % Stability behaviour can be determined at the last curve point
            else                                                % In all other cases (e.g. two BFPs are following each other directly):
                idx_stab = 'plot_in_dark_grey';                 % Plot in dark grey since stability behaviour cannot be determined securely
            end
        else                                                    % If two idx points are not following each other ...
            idx_stab = idx(k)+1;                                % ... stability is determined at idx(k)+1, since idx(k) can be a BFP
        end                                                     % (The BFP can be stable as well as unstable due to being calculated by interval switching)

        % Set the color depending on the stability behaviour
        if ischar(idx_stab)
            options.color = obj.plot_color.dark_grey;
        elseif n_unstable(idx_stab) == 0
            options.color = 'b';
        else
            options.color = 'r';
        end
           
        plot(mu(idx(k):idx(k+1)),s_out(idx(k):idx(k+1)),'color',obj.custom_color(options),'linestyle',linestyle); 

    end

end


% Plot the bifurcation labels
if any(bifurc_idx_local ~= 0)

    local_name = obj.bifurcation.bifurcation(bfp_bool);     % This gives the name of the bifurcation contained in the respective range of options.index.
    
    for k = 1:numel(bifurc_idx_local)
        text(mu(bifurc_idx_local(k)),s_out(bifurc_idx_local(k)),append('$\quad$',local_name{k}),'HorizontalAlignment','left','Fontsize',12,'Interpreter','latex');
        plot(mu(bifurc_idx_local(k)),s_out(bifurc_idx_local(k)),'Marker','.','Markersize',15,'MarkerFaceColor','k','MarkerEdgeColor','k');
    end

end


% Title, labels and x-limits
title(TtlStr,'Interpreter','latex');
xlabel('Continuation parameter $\mu$','Interpreter','latex');
ylabel(ystr,'Interpreter','latex');
xlim(DYN.opt_cont.mu_limit)



%% Output
varargout{1,1} = s_out;
varargout{1,2} = mu;


end