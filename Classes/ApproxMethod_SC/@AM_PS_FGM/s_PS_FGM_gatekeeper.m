%Gatekeeper function for the AM_PS_FGM subclass. In here, all input
%parameters are checked, before further processing. This static method is
%called from the static method s_AM_gatekeeper
%
%@GC:                   object of Gatekeeper class 
%@system:               user supplied option structure for the system. Needed
%here for checking the correct dimension of the initial value
%@opt_sol:              user supplied option structure for solution. Needed
%here for checking the phasecondition
%@opt_approx_method:       user supplied option structure for the solution
%method
%@opt_init_sol:         user supplied option structure for initializing a
%solution

function s_PS_FGM_gatekeeper(GC,system,opt_sol,opt_approx_method,opt_init)

    opt_approx_method_mandatory_fieldnames  = {};                                                                                  %mandatory fieldsnames in the options super structure
    opt_approx_method_allowed_fieldnames    = {'n_fft','phasecond','error_control','error_limit','ec_iter_max','n_hh_max'};         %allowed fieldsnames in the options super structure

    opt_init_mandatory_fieldnames_1  = {'c0','cmatrix','smatrix','hmatrix'};                      %mandatory fieldsnames in the options super structure
    opt_init_allowed_fieldnames_1    = {'c0','cmatrix','smatrix','hmatrix'};                      %allowed fieldsnames in the options super structure

    opt_init_mandatory_fieldnames_2   = {'fc0','hmatrix'};                          %mandatory fieldsnames in the options super structure
    opt_init_allowed_fieldnames_2     = {'fc0','hmatrix'};                          %allowed fieldsnames in the options super structure


    phasecond_allowed_fieldsvalues = {'poincare','int_poincare'};  %allowed fieldsnames for phasecondition values
    %% Check the opt_approx_method structure
  
    GC.check_fields(opt_approx_method,'opt_approx_method',opt_approx_method_mandatory_fieldnames,opt_approx_method_allowed_fieldnames);
    GC.check_fields(opt_init,'opt_init',opt_init_mandatory_fieldnames_1,opt_init_allowed_fieldnames_1,opt_init_mandatory_fieldnames_2,opt_init_allowed_fieldnames_2);
    GC.speak();

    %%%%%%%%%%%%%% opt_init CHECK %%%%%%%%%%%%%%%%%
    %% check the entries data types
    %Check the mandatory fields first (these are defintively present)
    %%%%%%%%%%%%%%%%%%%%
    %Get these values for the size check
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check the hmatrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_hh    = numel(opt_init.hmatrix);
    GC.check_data(opt_init.hmatrix,'opt_init.hmatrix','double','vector',[]);        %This must be a vector, since the hmatrix always contains a 0 and at least a 0 or a 1
 
    tmp = opt_init.hmatrix(:).';
    if ~all(tmp(:,1)==zeros(1,1))
        GC.error_msg{1,end+1} = append('The first frequency combination of  opt_init.hmatrix must be 0. However, it is ',num2str(tmp(1,1)));
    end 
    if ~any(tmp==1)
        GC.error_msg{1,end+1} = append('Your opt_init.hmatrix must at least contain the frequency 1.');
    end

    %Check that there are no doubled frequencies
    [~,IA,IC] = unique(abs(opt_init.hmatrix)); %abs ensures that no negative and positive doubles are contained
    
    if ~(numel(IA)==numel(IC))
       GC.error_msg{1,end+1} =  append('Your opt_init.hmatrix contains ',num2str(numel(IC)-numel(IA)),' doubled entry/entries, but must be unique.');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %This whole block is necessary, since the number of initial conditions can be equal OR smaller 
    %than the number of higher harmonics (the missing ICs will be guessed by CoSTAR)
    if isfield(opt_init,'c0')
        GC.check_data(opt_init.c0,'opt_init.c0','double',{'scalar','vector'},[]);
        GC.speak();

        %Check data type and that matrices are matrices
        GC.check_data(opt_init.cmatrix,'opt_init.cmatrix','double',{'vector','matrix'},[]);
        GC.check_data(opt_init.smatrix,'opt_init.smatrix','double',{'vector','matrix'},[]);        %This must be a vector, since the hmatrix always contains a 0 and at least a 0 or a 1
    
        %Check the correct size of the matrices: the number of entries must be equal or smaller than the number of harmonics (in the latter case Fourier coefficient values are guessed)
        if ~(size(opt_init.cmatrix,1)==system.dim)
            GC.error_msg{1,end+1} = append('Your input opt_init.cmatrix matrix has the dimension [',num2str(size(opt_init.cmatrix,1)), ...
                ',',num2str(size(opt_init.cmatrix,2)),'] but must be [',num2str(system.dim), ...
                ',n] with n <=',num2str(n_hh-1));
        end
        if (size(opt_init.cmatrix,2)>n_hh-1)
             GC.error_msg{1,end+1} = append('Your input opt_init.cmatrix matrix has the dimension [',num2str(size(opt_init.cmatrix,1)), ...
                ',',num2str(size(opt_init.cmatrix,2)),'] but must be [',num2str(system.dim), ...
                ',n] with n <=',num2str(n_hh-1),' with ',num2str(n_hh-1),' being the number of harmonics without constant term.');
        end

        if ~(size(opt_init.smatrix,1)==system.dim)
            GC.error_msg{1,end+1} = append('Your input opt_init.smatrix matrix has the dimension [',num2str(size(opt_init.smatrix,1)), ...
                ',',num2str(size(opt_init.smatrix,2)),'] but must be [',num2str(system.dim), ...
                ',n] with n <=',num2str(n_hh-1));
        end

        if (size(opt_init.smatrix,2)>n_hh-1)
             GC.error_msg{1,end+1} = append('Your input opt_init.smatrix matrix has the dimension [',num2str(size(opt_init.smatrix,1)), ...
                ',',num2str(size(opt_init.smatrix,2)),'] but must be [',num2str(system.dim), ...
                ',n] with n <=',num2str(n_hh-1),' with ',num2str(n_hh-1),' being the number of harmonics without constant term.');
        end

        if ~(size(opt_init.smatrix,2)==size(opt_init.cmatrix,2))
            GC.error_msg{1,end+1} = append('Your input opt_init.smatrix and opt_init.cmatrix must have the same size.');
        end

    else

        GC.check_data(opt_init.fc0,'opt_init.fc0','double','vector',[]);        %-1 since n_hh also contains the constant 0 frequency.
        if (numel(opt_init.fc0)>((2.*(n_hh-1)+1)*system.dim))||(numel(opt_init.fc0)<(3*system.dim))
                GC.error_msg{1,end+1} = append('Your input opt_init.fc0 vector has the dimension [',num2str(numel(opt_init.fc0)), ...
                ',1] but must be [n,1] with n <= ',num2str(((2.*n_hh-1)*system.dim)),' and n >= ',num2str((3*system.dim)),'.');
        end

        if ~(mod((numel(opt_init.fc0)-system.dim)./system.dim,2)==0)
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
       if max(opt_init.hmatrix) > opt_approx_method.n_fft/2
         
         GC.error_msg{1,end+1} = append('Your supplied value of opt_approx_method.n_fft = ', num2str(opt_approx_method.n_fft) ,' in combination with the maximum value of opt_init.hmatrix = ', num2str(max(opt_init.hmatrix)),' violates the Nyquist-Shannon Theorem. Consider raising opt_approx_method.n_fft.');   
       end
           GC.speak();
    end

%The s_DYN_gatekeeper already checked, that there is either non_auto_freq
%OR auto_freq
    if isfield(opt_approx_method,'phasecond')&&isfield(opt_sol,'non_auto_freq')
        GC.error_msg{1,end+1} = append('You are trying to solve a non-autonomous ODE, since you set the option opt_sol.non_auto_freq. In that case you are not allowed to set a phasecondition via opt_approx_method.phasecond.');   
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
    %Use the GC.check_data function here
    %% Do logical checks




   

end
