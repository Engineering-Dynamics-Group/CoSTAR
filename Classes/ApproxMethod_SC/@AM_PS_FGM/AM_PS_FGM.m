% Class AM_PS_FGM provides functions which return the residual for
% Fourier-Galerkin algorithms. AM_PS_FGM is a subclass of ApproxMethod
%
% @options:             struct contain all options for calculating solution datatype of all fields of @options have to be structs itself
%
%                       Mandatory fields (datatype: struct):
%                       'opt_init':           Mandatory fields of 'system':
%                                               --> 'hmatrix':                  Dimension of the state space

%                                               --> 'c0':                       0-th order Fourier coefficient. Size: [dim,1]
%                                               --> 'cmatrix':                  Cosine-forier coefficient matrix. Size [dim,numel(hmatrix)-1]
%                                               --> 'smatrix':                  Sine-forier coefficient matrix. Size [dim,numel(hmatrix)-1]
%                                               !!!!!OR!!!!!!
%                                               --> 'FC0':                      Fourier coefficient vector. Not recommended: Order of cosine and sine coeffs could be mixed up.


classdef AM_PS_FGM < ApproxMethod

    properties
        %Define here the properties, which are supplied in the
        %opt_approx_method 
        
        c0                          %Fourier coefficient vector or scalar for zeroth order higher harmonic
        cmatrix                     %Matrix of cosine Fourier coefficients - possibily provided by the user
        smatrix                     %Matrix of sine Fourier coefficients - possibily provided by the user
        fc0                         %Initial Fourier coefficient vector: Alternative to providing the matrices c0, smatrix, cmatrix

        hmatrix                     %higher harmonic vector
        
        %Define here the properties, which are supplied in the
        %opt_init 

        n_fft                       %Number of FFT points used for Fourier-Series build-up
        phasecond = 'poincare';     %Phase condition 
        

        %Everything for error control
        error_limit  = [1e-3,1e-1];     %Limit for the spectral error 
        ec_iter_max  = 10;              %Maximal iteration number for error control
        n_hh_max     = Inf;             %Maximum number of higher harmonics (for automatic increase by error_control)
        ec_prop_save = struct('iv',[],'hmatrix',[],'n_fft',[]);     %saves properties that are modified by the error control. If error control fails, the properties are resetted to the values stored in this struct
        
    end

        %Inherited Properties from Superclass 
        %res function_handle
        %iv initial value vector
        %error_control switch for the error control feature

    properties(Access = private)
        p_n_hh                                                      %number of higher harmonics
        p_arg_val                                                   %argument vector of base frequency (length is p_n_hh) for evaluation
        p_chf                                                       %complex harmonic function for building up the fourier series
        p_hmatrix_old                                               %hmatrix of the last iteration. This is needed for the error_control 


    end
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    methods(Static)                                                        %Static: Method can be called without creating an object of class AM_PS_FGM;

        s_PS_FGM_gatekeeper(GC,system, opt_sol,opt_approx_method,opt_init);%Gatekeeper method, which is called by the static s_AM_gatekeeper method  to check the inputs at the beginning 
        help_struct = s_help_opt_approx_method_PS_FGM();                     %Help struct for the opt_approx_method option structure 
        help_struct = s_help_opt_init_PS_FGM();                              %Help struct for the opt_init option structure

    end
    %%%%%%%%%%%%%%%
    methods
        %% Constructor
        function obj = AM_PS_FGM(DYN)

            obj.error_control = 'on';          %Switch for using the error_control feature - overwrite the default value from the superclass ApproxMethod. Needs to be set before updateoptions

            %Set important properties (POSITION BEFORE UPDATEOPTIONS IS CRITICAL)
            obj.p_arg_val = 0:(2*pi/2^6):(2*pi-2*pi/2^6);
            obj.hmatrix = [0,1];                            %Dummy property - gets replaced by updateoptions
            obj.n_fft = 2^6;                                %Needs to be set here after p_arg_val, since this triggers the set method.
            
            
            
            obj = updateoptions(obj,DYN.opt_approx_method);                   %assign all fields in opt_approx_method to class properties 
            obj = updateoptions(obj,DYN.opt_init);                         %assign all fields in opt_init to class properties

            obj.error_limit = sort(obj.error_limit,'ascend');          %assure that the error_limit array is ordered in increasing direction
       
            %Calculate a few necessary parameters
            obj.hmatrix = abs(obj.hmatrix(:).');                           %Assure that hmatrix has the right orientation. Only positive frequencies allowed

            obj = obj.getIV(DYN);                                          %Set initial value (Has to be set here, because residual accesses iv
     
            %Adapt n_FFT if not supplied, and if default value is not sufficient. (if it supplied by the user it was already checked in the gatekeeper)
            if isfield(DYN.opt_init,'n_fft')
                if max(obj.hmatrix) > obj.n_fft/2
                    tmp = nextpow2(2.*max(obj.hmatrix)+1);
                    obj.n_fft = 2^tmp;
                     warning(append('Your highest harmonic supplied in hmatrix violates the Nyquist-Shannon criterium in combination with the default value of n_fft = 2^6. The value was adapted accordingly to 2^',num2str(tmp),'. Consider supplying a sufficient value via opt_approx_method.n_fft.'));
                end
            end
        end

        %% Set Methods: These set methods guarantee that the value of p_n_hh, p_chf and p_arg_val are always updated/refreshed values.

        function set.hmatrix(obj,value)
            obj.hmatrix = value;
            obj.p_n_hh  = numel(value);
            obj.p_chf   = exp(1i.*value'*obj.p_arg_val);
            if max(obj.hmatrix,[],'all') > obj.n_fft/2-1
                    tmp = nextpow2(2.*max(obj.hmatrix)+1);
                    warning(append('Your highest harmonic in hmatrix violates the Nyquist-Shannon criterium in combination with the value of n_fft = 2^',num2str(log2(obj.n_fft)),'. The value was adapted accordingly to 2^',num2str(tmp),'.'));
                    obj.n_fft = 2^tmp;
            end
        
        end

        function set.n_fft(obj,value)
            obj.n_fft       = value;
            obj.p_arg_val   = 0:(2*pi/value):(2*pi-2*pi/value);
            obj.p_chf       = exp(1i.*obj.hmatrix'*obj.p_arg_val);
        end
       
        %% Interface Methods
        %Interface methods: This is an abstract method and must be defined
        %Use this interface to pass information between the continuer
        %algorithm and your AM subclass
        %@CON: Continuation class object
        function obj = IF_up_res_data(obj,CON)                                     
        
            obj.iv = CON.yp(1:(end-1));             %update the current initial condition. Used for the poincare phase condition.

        end
       
        %Function for estimating the spectral error
        %@y1:   solution curve point vector
        %@DYN:  DynamicalSystem class object
        %@err:  spectral error
        function err = IF_estimate_error(obj,y1,DYN)    

                proj = obj.residuum_projection(y1,DYN);
                err = sum(vecnorm(abs(proj),2,1))./vecnorm(y1,2);           %Normalised Norm of errror in frequency space 

        end

        [yp,flag] = IF_increase_discretization(obj,CON,DYN);                %Method, which increases the number of harmonics by 1 
        [yp,flag] = IF_decrease_discretization(obj,CON,DYN);                %Method, which decreasses the number of harmonics by 1
        varargout = IF_update_sol_dim(obj,DYN,new_dim, varargin);           %Method which updates previous curve point varargin/out{1,1} or tangent/secant varargin/out{1,2} properties used in according to a new solution space dimension new_dim
        IC        = getIC(obj,y,DYN,n_shoot);                                       %Method for getting an initial point in state space on the periodic solution orbit             

        %% Functions for residuum and computation
        res = PS_FGM_residuum(obj,y,DYN);
        obj = getIV(obj,DYN);
        
    end

    methods(Access = private)

        ph   = phase_condition(obj,FCtemp,DYN);                            %Method for the phase_condition
        proj = residuum_projection(obj,y,DYN);                             %Method for projecting the residuum of the system into the frequency space
        [s,hmatrix] = sort_guess_FC(obj,DYN,FC0);                          %sorting and possibly guessing new Fourier coefficients according to the higher harmonics

    end
end