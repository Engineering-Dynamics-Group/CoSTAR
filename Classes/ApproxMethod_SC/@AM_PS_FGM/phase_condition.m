%This function is a method of subclass AM_PS_FGM. This function computes
%a phase condition for the solution of autonomus periodic systems
%
%@obj:      ApproxMethod subclass object
%@FCtemp:   Current Fourier coefficient vector 
%@DYN:      DynamicalSystem object 
%
%@ph:       phase condition/ equation, which gets added to the residuum
%equation.


function ph = phase_condition(obj,FCtemp,DYN)


switch obj.phasecond

    case 'poincare' %This one works at the moment.

        x0 = obj.iv;
        dim = DYN.dim;
        Gamma0 = obj.iv(end);
        p_n_hh = obj.p_n_hh;

        FC0 = [x0(1:dim);x0((dim+1):(p_n_hh)*dim)-1i.*x0(((p_n_hh)*dim+1):(end-1))];    %Assemble complex Fourrier vector
        FCtemp0 = reshape(FC0,[dim,p_n_hh]);

        FS0 = real(FCtemp0*obj.p_chf);
        FS = real(FCtemp*obj.p_chf);
        dFS0 = real(1i.*Gamma0.*FCtemp0*(obj.hmatrix.'.*obj.p_chf));

        ph = dFS0(:,1).'*(FS(:,1)-FS0(:,1));

    case 'int_poincare' %int_poincare evaluated in time domain: Integral is computed analytically 

        x0 = obj.iv;
        dim = DYN.dim;
        p_n_hh = obj.p_n_hh;

        %!!! This needs attention w.r.t. matrix size, if adaption of higher harmonics is incorporated!!!
        FC0 = [x0((dim+1):(p_n_hh)*dim)-1i.*x0(((p_n_hh)*dim+1):(end-1))];          %Assemble complex Fourrier vector of the previously known solution (or start point of Newton-Scheme at the first curve point)
        FCtemp0 = reshape(FC0,[dim,p_n_hh-1]);                                      %Reshape to matrix
        ph = -sum(obj.hmatrix(2:end).*sum(imag(conj(FCtemp(:,2:end)).*FCtemp0),1)); %Analytical solution of the integral Poincare condition (possible for FGM)

end

end