%Gatekeeper function for the AM_template subclass. In here, all input
%parameters are checked, before further processing. This static method is
%called from the static method s_AM_gatekeeper
%
%@GC:                   object of Gatekeeper class 
%@opt_approx_method:       user supplied option structure for the solution
%method
%@opt_init_sol:         user supplied option structure for initializing a
%solution

function s_QPS_FGM_gatekeeper(GC,system,opt_sol,opt_approx_method,opt_init)

    opt_approx_method_mandatory_fieldnames  = {};                                                                              %mandatory fieldsnames in the options super structure
    opt_approx_method_allowed_fieldnames    = {'n_fft','phasecond','error_control','error_limit','ec_iter_max','n_hh_max'};    %allowed fieldsnames in the options super structure

    opt_init_mandatory_fieldnames_1  = {'c0','cmatrix','smatrix','hmatrix'};                                                %mandatory fieldsnames in the options super structure
    opt_init_allowed_fieldnames_1    = {'c0','cmatrix','smatrix','hmatrix'};                                                %allowed fieldsnames in the options super structure

    opt_init_mandatory_fieldnames_2   = {'fc0','hmatrix'};                                      %mandatory fieldsnames in the options super structure
    opt_init_allowed_fieldnames_2     = {'fc0','hmatrix'};                      %allowed fieldsnames in the options super structure


    phasecond_allowed_fieldsvalues = {'poincare','int_poincare'};                            %allowed fieldsnames for phasecondition values
    %% Check the opt_approx_method structure
  
    GC.check_fields(opt_approx_method,'opt_approx_method',opt_approx_method_mandatory_fieldnames,opt_approx_method_allowed_fieldnames);
    GC.check_fields(opt_init,'opt_init',opt_init_mandatory_fieldnames_1,opt_init_allowed_fieldnames_1,opt_init_mandatory_fieldnames_2,opt_init_allowed_fieldnames_2);
    GC.speak();

    %%%%%%%%%%%%%% opt_init CHECK %%%%%%%%%%%%%%%%%
    %% check the entries data types
    %Check the mandatory fields first (these are defintively present)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Check the hmatrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [n_hh,idx] = sort(size(opt_init.hmatrix));                                    %The first entry of n_hh is the smaller dimension must always be 2! The second determines the number of higher
    
    %harmonics  must at least be 3
    %Check hmatrix is a matrix
    GC.check_data(opt_init.hmatrix,'opt_init.hmatrix','double','matrix',[]);        %This must be a matrix, since the hmatrix always contains a 0,0 a 0,1 and a 1,0 
    GC.speak();
    %Check its dimensions
    if ~(n_hh(1)==2)
        GC.error_msg{1,end+1} = append('Only a 2-D quasi-periodic method is implemented at the moment: The smallest dimension of opt_init.hmatrix must be 2, but is at the moment ',num2str(n_hh(1)));
    end
    if  n_hh(2)<3
        GC.error_msg{1,end+1} = append('Your opt_init.hmatrix be a 2x"number_of_harmonics" sized matrix. At least 3 harmonics (including constant term) must be present. However, your largest dimension is ',num2str(n_hh(2)));
    end
    GC.speak();
    
    %Check the first three entries:
    tmp = opt_init.hmatrix;
    tmp = permute(tmp,idx);     %Make sure tmp is a [2 x numb_hh]-sized matrix
    
    if ~all(tmp(:,1)==zeros(2,1))
        GC.error_msg{1,end+1} = append('The first frequency combination of  opt_init.hmatrix must be [0;0]. However, it is [',num2str(tmp(1,1)),',',num2str(tmp(2,1)),']');
    end 
    if ~any(prod((tmp==[1;0]),1))
        GC.error_msg{1,end+1} = append('Your opt_init.hmatrix must at least contain the frequency combination [1;0].');
    end
    if ~any(prod((tmp==[0;1]),1))
        GC.error_msg{1,end+1} = append('Your opt_init.hmatrix must at least contain the frequency combination [0;1].');
    end
    GC.speak();

    %Check that there are no doubled frequencies
    [~,IA,IC] = unique(permute(opt_init.hmatrix,idx)','rows'); %permute + transpose ensures, that hmatrix is alway n_hh x 2 matrix (better for this case)
    
    if ~(numel(IA)==numel(IC))
       GC.error_msg{1,end+1} =  append('Your opt_init.hmatrix contains ',num2str(numel(IC)-numel(IA)),' doubled entry/entries, but must be unique.');
    end

    GC.speak();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check the initial values for the Fourier coefficients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This whole block is necessary, since the number of initial conditions can be equal OR smaller 
    %than the number of higher harmonics (the missing ICs will be guessed by CoSTAR)
    
    if isfield(opt_init,'c0')
        GC.check_data(opt_init.c0,'opt_init.c0','double',{'scalar','vector'},[]);
        GC.speak();

        %Check data type and that matrices are matrices
        GC.check_data(opt_init.cmatrix,'opt_init.cmatrix','double',[],[]);
        GC.check_data(opt_init.smatrix,'opt_init.smatrix','double',[],[]);        %This must be a matrix since the hmatrix always contains a [0,0], a [1,0] and a [0,1]
    
        %Check the correct size of the matrices: the number of entries must be equal or smaller than the number of harmonics (in the latter case Fourier coefficient values are guessed)
        if ~(size(opt_init.cmatrix,1)==system.dim)
            GC.error_msg{1,end+1} = append('Your input opt_init.cmatrix matrix has the dimension [',num2str(size(opt_init.cmatrix,1)), ...
                ',',num2str(size(opt_init.cmatrix,2)),'] but must be [',num2str(system.dim), ...
                ',n] with n <=',num2str(n_hh(1,2)-1));
        end
        if (size(opt_init.cmatrix,2)>n_hh(1,2)-1)
             GC.error_msg{1,end+1} = append('Your input opt_init.cmatrix matrix has the dimension [',num2str(size(opt_init.cmatrix,1)), ...
                ',',num2str(size(opt_init.cmatrix,2)),'] but must be [',num2str(system.dim), ...
                ',n] with n <=',num2str(n_hh(1,2)-1),' with ',num2str(n_hh(1,2)-1),' being the number of harmonics without constant term.');
        end

        if ~(size(opt_init.smatrix,1)==system.dim)
            GC.error_msg{1,end+1} = append('Your input opt_init.smatrix matrix has the dimension [',num2str(size(opt_init.smatrix,1)), ...
                ',',num2str(size(opt_init.smatrix,2)),'] but must be [',num2str(system.dim), ...
                ',n] with n <=',num2str(n_hh(1,2)-1));
        end

        if (size(opt_init.smatrix,2)>n_hh(1,2)-1)
             GC.error_msg{1,end+1} = append('Your input opt_init.smatrix matrix has the dimension [',num2str(size(opt_init.smatrix,1)), ...
                ',',num2str(size(opt_init.smatrix,2)),'] but must be [',num2str(system.dim), ...
                ',n] with n <=',num2str(n_hh(1,2)-1),' with ',num2str(n_hh(1,2)-1),' being the number of harmonics without constant term.');
        end

        if ~(size(opt_init.smatrix,2)==size(opt_init.cmatrix,2))
            GC.error_msg{1,end+1} = append('Your input opt_init.smatrix and opt_init.cmatrix must have the same size.');
        end

    else

        GC.check_data(opt_init.fc0,'opt_init.fc0','double','vector',[]);        %-1 since n_hh also contains the constant 0 frequency.
        if (numel(opt_init.fc0)>((2.*(n_hh(1,2)-1)+1)*system.dim))||(numel(opt_init.fc0)<(5*system.dim))
                GC.error_msg{1,end+1} = append('Your input opt_init.fc0 vector has the dimension [',num2str(numel(opt_init.fc0)), ...
                ',1] but must be [n,1] with n <= ',num2str(((2.*n_hh(1,2)+1)*system.dim)),' and n >= ',num2str((5*system.dim)),'.');
        end

        if ~(mod((numel(opt_init.fc0)-system.dim)./(system.dim),2)==0)
            GC.error_msg{1,end+1} = append('Your input opt_init.fc0 vector must have an equal number of initial conditions for the sine and cosine terms.');
        end 

    end

    GC.speak();

  
    %%%%%%%%%%%%%% opt_approx_method CHECK %%%%%%%%%%%%%%%%%
    %Check the optional fields now 
    %%%%%%%%%%%%%%%%%%%%
    if isfield(opt_approx_method,'n_fft')
       GC.check_data(opt_approx_method.n_fft,'opt_approx_method.n_FFT','double','scalar',[]); 
       GC.speak();
       %Check if n_FFT is in the power of two:
       if ~(mod(log(opt_approx_method.n_fft)/log(2),1)==0)
         GC.error_msg{1,end+1} = append('There is a problem with your value of opt_approx_method.n_fft: The current value is ', num2str(opt_approx_method.n_fft), ' but must be a power of two (e.g. 2^6)');   
       end
       if max(opt_init.hmatrix,[],'all') > opt_approx_method.n_fft/2
         
         GC.error_msg{1,end+1} = append('Your supplied value of opt_approx_method.n_fft = ', num2str(opt_approx_method.n_fft) ,' in combination with the maximum harmonic of opt_init.hmatrix = ', num2str(max(opt_init.hmatrix,[],'all')),' violates the Nyquist-Shannon Theorem. Consider raising opt_approx_method.n_fft.');   
       end
           GC.speak();
    end
    
    if isfield(opt_approx_method,'phasecond')&&isfield(opt_sol,'non_auto_freq')
        if ~isfield(opt_sol,'auto_freq')
            GC.error_msg{1,end+1} = append('You are trying to solve a non-autonomous ODE, since you set the option opt_sol.non_auto_freq. In that case you are not allowed to set a phasecondition via opt_approx_method.phasecond.');   
        end
    end
    if isfield(opt_approx_method,'phasecond')&&isfield(opt_sol,'auto_freq'); GC.check_data(opt_approx_method.phasecond, 'opt_approx_method.phasecond','char',[],phasecond_allowed_fieldsvalues); end
  
    if isfield(opt_approx_method,'error_control'); GC.check_data(opt_approx_method.error_control, 'opt_approx_method.error_control','char',[],{'on','off'}); end
  
    if isfield(opt_approx_method,'error_limit'); GC.check_data(opt_approx_method.error_limit, 'opt_approx_method.error_limit','double',[1,2],'positive'); end
    
    GC.speak();


     if isfield(opt_approx_method,'ec_iter_max') 
         GC.check_data(opt_approx_method.ec_iter_max, 'opt_approx_method.ec_iter_max','double','scalar','positive'); 
         GC.speak();
         if ~(isreal(opt_approx_method.ec_iter_max) && rem(opt_approx_method.ec_iter_max,1)==0)
            GC.error_msg{1,end+1} = append('Your parameter opt_approx_method.ec_iter_max = ',GC.my_join(opt_approx_method.ec_iter_max),' can only be a real integer.'); 
         end
     end

     if isfield(opt_approx_method,'n_hh_max') 
         GC.check_data(opt_approx_method.n_hh_max, 'opt_approx_method.n_hh_max','double','scalar','positive'); 
         GC.speak();
         if ~(isreal(opt_approx_method.n_hh_max) && rem(opt_approx_method.n_hh_max,1)==0)
            GC.error_msg{1,end+1} = append('Your parameter opt_approx_method.n_hh_max = ',GC.my_join(opt_approx_method.n_hh_max),' can only be a real integer.'); 
         end
     end

    GC.speak();
    %% Do logical checks

end
