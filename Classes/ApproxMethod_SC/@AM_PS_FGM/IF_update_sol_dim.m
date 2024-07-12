%This method updates the dimension of the solution data of past curve points, if the solution space dimension has been changed 
% due to the error_control method or within the iteration process for the bifurcation point (Stability class )
%
%@obj:              AM_PS_FGM class object
%@DYN:              Dynamical System class object
%@new_dim:          New dimension of the curve point. y0 and dy0 get updated to this dimension
%@varargin{1,1}:    y0 - solution curve point to be update in its dimension
%@varargin{1,2}:    dy0 - solution curve point tangent/secant to be update in its dimension

%@varargout{1,1}:    y0 -  updated solution curve point
%@varargout{1,2}:    dy0 - updated solution curve point tangent/secant 

%@CON: Continuation class object


function varargout = IF_update_sol_dim(obj,DYN,new_dim, varargin) 
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% y0 %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y0 = varargin{1,1};

    tmp = reshape(y0(1:end-1-DYN.n_auto),DYN.dim,[]);
    c0 = tmp(:,1);
    cmatrix = tmp(:,2:0.5*(numel(tmp)/DYN.dim+1));          %Gatekeeper has ensured that this computation of dimension is always possible
    smatrix = tmp(:,0.5*(numel(tmp)/DYN.dim+1)+1:end);

   
    if new_dim >numel(y0)
        
        numb = (new_dim-numel(y0))./2;
        varargout{1,1} = [c0(:);cmatrix(:);zeros(numb,1);smatrix(:);zeros(numb,1);y0(end-DYN.n_auto:end)];
   
    elseif new_dim<numel(y0)

        %p_hmatrix_old is the old hmatrix of the last iteration in the CON.error_control loop. You can't use S.hmatrix{1,end} here, since this is the state before
        %the error_control loop began. However, for multiple loops within the error_control, the state of the last loop of hmatrix is needed.
        hhmold = obj.p_hmatrix_old;

        tmp_hh = setdiff(hhmold,obj.hmatrix);
        if ~isempty(tmp_hh); idx = find(hhmold == tmp_hh); idx = idx-1;  end
        if isempty(idx); idx = size(cmatrix,2)-1; obj.hmatrix = hhmold(1,1:end-1); end %If something goes wrong - simply delete the last harmonic
        idx2= [1:(idx-1),(idx+1):numel(hhmold)-1];
        

        varargout{1,1} = [c0(:);reshape(cmatrix(:,idx2),[],1);reshape(smatrix(:,idx2),[],1);y0(end-DYN.n_auto:end)];
        
    else %This case occurs, if no new dimension could be found by the increase or decrease methods
        varargout{1,1} = y0;

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% dy0 %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if size(varargin,2) == 2 %Update of tangent is also necessary
        
        dy0 = varargin{1,2};

        tmp = reshape(dy0(1:end-1-DYN.n_auto),DYN.dim,[]);
        dc0 = tmp(:,1);
        dcmatrix = tmp(:,2:0.5*(numel(tmp)/DYN.dim+1));          %Gatekeeper has ensured that this computation of dimension is always possible
        dsmatrix = tmp(:,0.5*(numel(tmp)/DYN.dim+1)+1:end);

        if new_dim>numel(dy0)

            varargout{1,2} = [dc0(:);dcmatrix(:);zeros(numb,1);dsmatrix(:);zeros(numb,1);dy0(end-DYN.n_auto:end)];

        elseif new_dim<numel(dy0)

            varargout{1,2} = [dc0(:);reshape(dcmatrix(:,idx2),[],1);reshape(dsmatrix(:,idx2),[],1);dy0(end-DYN.n_auto:end)];

        else %This case occurs, if no new dimension could be found by the increase or decrease methods
            varargout{1,2} = dy0;

        end
    end



end





