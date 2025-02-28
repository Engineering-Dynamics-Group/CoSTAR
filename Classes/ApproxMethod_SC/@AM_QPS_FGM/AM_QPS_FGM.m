% Class AM_QPS_FGM provides functions which return the residual for the fourier-galerkin algorithm
% algorithms in the quasi-periodic case. AM_QPS_FGM is a subclass of ApproxMethod
%
classdef AM_QPS_FGM < ApproxMethod

    properties

        c0                              %Fourier coefficient vector or scalar for zeroth order higher harmonic
        cmatrix                         %Matrix of cosine Fourier coefficients - possibily provided by the user
        smatrix                         %Matrix of sine Fourier coefficients - possibily provided by the user
        fc0                             %Initial Fourier coefficient vector: Alternative to providing the matrices c0, smatrix, cmatrix

        hmatrix                         %higher harmonic vector
        
        n_fft = 2^6;                    %Number of FFT points used for Fourier-Series build-up
        phasecond = 'int_poincare'      %Phase condition 

        %Everything for error control
        error_limit = [1e-3,1e-1];      %Limit for the spectral error 
        ec_iter_max = 10;               %Maximal iteration number for error control
        n_hh_max      = Inf;          %Maximum number of higher harmonics (for automatic increase by error_control)
    end
    
    properties(Access=private)
        p_n_hh                      %Number of higher harmonics
        p_arg_val                   %Argument vector of base frequency (length is p_n_hh) for evaluation
        p_chf                       %complex harmonic function for building up the fourier series
        p_hmatrix_old               %hmatrix of the last iteration. This is needed for the error_control 


    end
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    methods(Static)                                                        %Static: Method can be called without creating an object of class Shoot;

        s_QPS_FGM_gatekeeper(GC,system,opt_sol,opt_approx_method,opt_init);         %Gatekeeper method, which is called by the static s_AM_gatekeeper method  to check the inputs at the beginning 
        help_text = s_help_opt_approx_method_QPS_FGM();                             %Help text for the opt_approx_method option structure 
        help_text = s_help_opt_init_QPS_FGM();                                      %Help text for the opt_init option structure 

    end
    %%%%%%%%%%%%%%%
    methods
        %% Constructor
        function obj = AM_QPS_FGM(DYN)  

            obj.error_control = 'on';          %Switch for using the error_control feature - overwrite the default value from the superclass ApproxMethod. Needs to be set before updateoptions

            %Set important properties (POSITION BEFORE UPDATEOPTIONS IS CRITICAL)
            tmp  = 0:(2*pi/2^6):(2*pi-2*pi/2^6);                                %evaluation theta_values
            obj.p_arg_val  = [reshape(repmat(tmp' ,  [1  ,2^6]),[1,(2^6)^2]);                       
                               reshape(repmat(tmp   ,[2^6,1  ]),[1,(2^6)^2])]; 

            obj.hmatrix = [0,1,0;0,0,1]; %Dummy property - gets replaced by updateoptions
            obj.n_fft = 2^6;    %Needs to be set here after p_arg_val, since this triggers the set method.


            obj = updateoptions(obj,DYN.opt_approx_method);                   %updateoptions method is inherited from ApproxMethod
            obj = updateoptions(obj,DYN.opt_init);                          %updateoptions method is inherited from ApproxMethod
       
            %% Check the higher harmonics matrix to be ok
            %Calculate a few necessary parameters
            [~,idx] = sort(size(obj.hmatrix));                             %These two command lines assure, that the first dimension of the hmatrix is the smaller one
            obj.hmatrix = permute(obj.hmatrix,idx);                        %(Gatekeeper checked that the correct dimension are present).
            obj.p_n_hh = size(obj.hmatrix,2);                              %Number of higher harmonics (0 for constant term is included)
            
            %Check that the second row only contains positive integers (as per Definition). If not-swap the sign.
            idx = find(obj.hmatrix(2,:)<0);
            obj.hmatrix(:,idx) = -obj.hmatrix(:,idx);
                
            obj = obj.getIV(DYN);                                          %Set initial value (Has to be set here, because residual accesses iv)
                                                                           %getIV also sorts the hmatrix according to ascending harmonics and the initial values accordingly
     
            %Adapt n_FFT if not supplied, and if default value is not sufficient. (if it supplied by the user it was already checked in the gatekeeper)
            if isfield(DYN.opt_init,'n_fft')
                if max(obj.hmatrix,[],'all') > obj.n_fft/2-1
                    tmp = nextpow2(2.*max(obj.hmatrix)+1);
                    obj.n_fft = 2^tmp;
                    warning(append('Your highest harmonic supplied in hmatrix violates the Nyquist-Shannon criterium in combination with the default value of n_fft = 2^6. The value was adapted accordingly to 2^',num2str(tmp),'. Consider supplying a sufficient value via opt_approx_method.n_fft.'));
                end
            end
 
        end
        
         %% Set Methods: These set methods guarantee that the value of p_n_hh, p_chf and p_arg_val are always updated/refreshed values.
        function obj = set.hmatrix(obj,value)
            obj.hmatrix = value;
            obj.p_n_hh  = size(value,2);
            %How does this work?: I am evaluating here the arguments in the cosine or sine functions: cos(H1*theta1+H2*theta2). hmatrix'*p_freq_val is the scalar product
            %The definition of p_freq_val (similar to the meshgrid function) ensures, that p_chf is evaluated at every point of the discretized torus. 
            obj.p_chf        = exp(1i.*value'*obj.p_arg_val);        %complex harmonic functions: e.g. real(obj.p_chf(2,:)) gives the cosine with the first higher harmonic defined in obj.hmatrix (first one is 0)
                                                                            %imag(obj.p_chf(2,:)) correspondingly gives the same for the sine function
            if max(obj.hmatrix,[],'all') > obj.n_fft/2-1
                    tmp = nextpow2(2.*max(obj.hmatrix,[],'all')+1);
                    warning(append('Your highest harmonic in hmatrix violates the Nyquist-Shannon criterium in combination with the value of n_fft = 2^',num2str(log2(obj.n_fft)),'. The value was adapted accordingly to 2^',num2str(tmp),'.'));
                    obj.n_fft = 2^tmp;
            end
        
        end

        function obj = set.n_fft(obj,value)
            obj.n_fft       = value;
            tmp  = 0:(2*pi/value):(2*pi-2*pi/value);  
            obj.p_arg_val   =  [reshape(repmat(tmp' ,  [1    ,value]),[1,(value)^2]);                       
                               reshape(repmat(tmp   ,  [value,1    ]),[1,(value)^2])]; 

            %How does this work?: I am evaluating here the arguments in the cosine or sine functions: cos(H1*theta1+H2*theta2). hmatrix'*p_freq_val is the scalar product
            %The definition of p_freq_val (similar to the meshgrid function) ensures, that p_chf is evaluated at every point of the discretized torus. 

            obj.p_chf        = exp(1i.*obj.hmatrix'*obj.p_arg_val);        %complex harmonic functions: e.g. real(obj.p_chf(2,:)) gives the cosine with the first higher harmonic defined in obj.hmatrix (first one is 0)
                                                                            %imag(obj.p_chf(2,:)) correspondingly gives the same for the sine function
        end
       
    


        %% Methods
        %Interface methods: This is an abstract method and must be defined
        %Use this interface to pass information between the continuer
        %algorithm and your AM subclass. 
        % @CON: Continuation class object
        function obj = IF_up_res_data(obj,CON)                                     
        
            %Exemplary usage:
            obj.iv = CON.yp(1:(end-1));                                         %update the current initial condition. Used for the poincare phase condition.

        end
    
        %Function estimates the current (spectral) error.
        %@y1:   current solution point vector
        %@DYN:  DynamicalSystem object
        %@err:  computed error
        function err = IF_estimate_error(obj,y1,DYN)    %function for estimatin the spectral error

                proj = obj.residuum_projection(y1,DYN);
                err =  sum(sum(vecnorm(abs(proj),2,3)))./vecnorm(y1,2);
           
          end

        [yp,flag] = IF_increase_discretization(obj,CON,DYN);                      %Method, which increases the number of harmonics by 1 
        [yp,flag] = IF_decrease_discretization(obj,CON,DYN);                      %Method, which decreasses the number of harmonics by 1
        varargout = IF_update_sol_dim(obj,DYN,new_dim, varargin)                  %Method which updates previous curve point varargin/out{1,1} or tangent/secant varargin/out{1,2} properties used in according to a new solution space dimension new_dim

       
        %% Functions for residuum and computation
        res = QPS_FGM_residuum(obj,y,DYN);
        obj = getIV(obj,DYN);  
        IC = getIC(obj,y,DYN,n_char_st);
    end

        methods(Access = private)

            ph  = phase_condition(obj,FCtemp,DYN);                           %Method for the phase_condition
            proj = residuum_projection(obj,y,DYN);                             %Method for projecting the residuum of the system into the frequency space
            [s,hmatrix] = sort_guess_FC(obj,DYN,FC0);                          %sorting and possibly guessing new Fourier coefficients according to the higher harmonics

        end


end