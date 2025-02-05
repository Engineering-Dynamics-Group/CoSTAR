% This static function creates an according solution object for the
% specific solution type and solution method
%
%@S:        Solution object
%@DYN:      DynamicalSystem object
%@AM:       ApproxMethod object

function ST = s_stability_selection(DYN,AM)


    %%%%%%%%%%%%% EQUILIBRIUM %%%%%%%%%%%%%
    if(strcmpi(DYN.sol_type,'equilibrium'))
    
        ST = ST_EQ(DYN);
        ST.calc_stability = @(y,J) ST.EQ_calc_stability(y,J,DYN);
        ST.test_functions = @(multipliers) ST.EQ_test_functions(multipliers,DYN);



    %%%%%%%%%%%%%%% PERIODIC %%%%%%%%%%%%%%
    elseif(strcmpi(DYN.sol_type,'periodic'))


        %%%%%%%%%%% SHOOTING %%%%%%%%%%%
        ST = ST_PS_SHM(DYN);                                                        % Stability computation of periodic solutions via SHM (possible for FDM and FGM as well!)
        ST.calc_stability = @(y,J) ST.PS_SHM_calc_stability(y,J,DYN,AM);            % Stability function
        ST.test_functions = @(multipliers) ST.PS_test_functions(multipliers,DYN);   % Method to evaluate the test functions


        %%%%%%%%%%%%% HILL %%%%%%%%%%%%% - to be implemented
        % ST.test_functions = @(multipliers) ST_PS_SHM.periodic_test_functions(multipliers,DYN);  % This is an abstract method and can be called without creating an object of ST_PS_SHM



    %%%%%%%%%%%% QUASI-PERIODIC %%%%%%%%%%%
    elseif(strcmpi(DYN.sol_type,'quasiperiodic'))
       

        %%%%%%%%%%% SHOOTING %%%%%%%%%%%
        if strcmpi(DYN.approx_method,'shooting')

            ST = ST_QPS_SHM(DYN);
            ST.calc_stability = @(y,J) ST.QPS_SHM_calc_stability(y,J,DYN,AM);           % Stability method for quasi-periodic shooting
            ST.test_functions = @(multipliers) NaN;                                     % There is no test function in the quasi-periodic case
        

        %%%%%%%%%%% FGM & FDM %%%%%%%%%%
        else

            ST = struct;        % Stability computation of quasi-periodic solutions only possible for SHM at the moment

        end

        
    end


end