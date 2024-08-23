% This function interpolates time integration data from an ODE Solver and makes a FFT
%
% @T: Timesteps from ODE solver
% @Z: Integration data from ODE solver
% @varargin: optional strings indicating window function (used to minimize leckage) and format of return
%             - Can be empty. Default window is 'rectangular' and default output format is 'polar'
%             - Available strings for output:
%                * 'polar':     f, |S| and angle(S)
%                * 'cartesian': f, real(S) and imag(S)
%             - Available windows:
%                * rectangular:   No change in input leads for non-periodic signals to leckage
%                * hann:          von-Hann window is quite universal. Good frequency resolution with moderate side-lobes
%                * humming:       Similar to the von-Hann but smaller first side-lobe. The remaining lobes are bigger
%                * blackmann:     Three term version. Wide main peak, but very good compression of side-lobes
%                * kaiser-bessel: Tries a compromise between accuary and compression. Similar to the blackmann window: for the same main ...
%                                 lobe width, higher first side lobes but smaller remaining lobes. Often reveals signals covered by noise.
%
% Example: costarFFT(T,Z,'window',name,'output',name);
%          costarFFT(T,Z,'window','hann','output','cartesian');
%
%
% Further explanations:
%
% Factor of 2 in S(2:end) = 2*S(2:end); (line 168):     CHECK LINES LATER ON
% The integral  F(y) = \int_(-\infty)^\infty f(x)exp(-2 \pi i x y)dx  can be
% rewritten as  F(y) = 2 \int_(0)^\infty f(x)exp(-2 \pi i x y)dx ,    which is the fourier signal we want, since we are only interested ...
% in the single sided spectrum. S(1) does not need to be multiplied by 2 because it is the amplitude of the constant part of the signal ...
% (frequency = 0), which does not have a complex conjugate pair.
%
% The "-" at -imag(S) (lines 172 and 176):       CHECK LINES LATER ON
% The Fourier Transform is defined as
% F(y) = \int_(-\infty)^\infty f(x)exp(-2 \pi i x y)dx = \int_(-\infty)^\infty f(x)[cos(2\pi x y)-i sin(2 \pi x y)]dx .
% The "-" in the exponent of the exponential function makes it necessary to multiply the result with a minus to get the correct sine amplitude.
%
%
% Copyright (c): Simon Baeuerle - University of Kassel - baeuerle@uni-kassel.de
% 
% V2: 02-2022: Bug-Fix: Z can now be a matrix and is correctly fast-fourier transformed.


function [f,varargout] = costarFFT(T,Z,varargin)

error_msg = ['Options must either be empty (default, output format polar, window rectangular) or contain the information of "window" (return '...
             'output polar) or and the information of "output" (default window: rectangular). Syntax is costarFFT(T,Z,"window",name,"output",name)'];
 

%% Initial Checks and presettings

%Convert all input arguments to lower case
for k=1:size(varargin,2)
    varargin{1,k} = lower(varargin{1,k});
end

%Check if inputs are non-empty
if isempty(T)||isempty(Z)
    error('Input variable T or Z are empty');
end


%% Set variables for output format and window  

if isempty(varargin)
    output = 'polar';
    window = 'rectangular';
    
elseif length(varargin) == 2
    
    if strcmpi(varargin{1,1},'output')
        output = varargin{1,2};
        window = 'rectangular';
    elseif strcmpi(varargin{1,1},'window')
        output = 'polar';
        window = varargin{1,2};
    else
        error(error_msg);
    end

elseif length(varargin) == 4
    
    if strcmpi(varargin{1,1},'output')
        output = varargin{1,2};
        if strcmpi(varargin{1,3},'window')
            window = varargin{1,4};
        else
            error(error_msg);
        end

    elseif strcmpi(varargin{1,1},'window')
        window = varargin{1,2}; 
        if strcmpi(varargin{1,3},'output')
            output = varargin{1,4};
        else
            error(error_msg);
        end
    else
        error(error_msg);
    end 

else
    error(error_msg);    
end


%% Interpolation of ODE data (in case not equally spaced)

[~,idx] = sort(size(Z),'descend');
Z = permute(Z,idx);                         % Make sure that Z is a column vector (expected by interp1 and FFT)

% t_step = (max(T)-min(T))/(length(T)-1);   % Time step of input signal -> moved down since interpolation produces a new time step when length(T) is uneven
% Ti = 0:t_step:(max(T)-min(T));            % Replaced by linspace to avoid numerical error
% 
% Zi = interp1((T-min(T)),Z, Ti);           % Interpolation is moved down to take into account that length(Z) can be uneven number
% 
% [~,idx] = sort(size(Zi),'descend');
% Zi = permute(Zi,idx);                     % Make sure Zi is a column vector (excpected by FFT) -> Not clear why again, thus commented
% 
% if rem(length(Zi),2) == 1                 % Make sure Zi has an even number of elements
%     Zi = Zi(1:(end-1),:);                 % Replaced by interpolation as removing one element can reduce the accuracy of the FFT
% end

if rem(length(T),2) == 1
    Ti = linspace(0,max(T)-min(T),length(T)-1)';    % Make sure Zi has even number of elements by interpolating on even number of query points
else
    Ti = linspace(0,max(T)-min(T),length(T))';      % T already has even number of elements
end
t_step = (max(Ti)-min(Ti))/(length(Ti)-1);          % Time step of input signal for FFT

Zi = interp1(T-min(T), Z, Ti, 'spline');            % Use spline interpolation. Every column of Z is treated separately


%% Apply window

switch lower(window)

    case 'rectangular'
        %Nothing needs to be done here
        
    case 'hann'
        N = length(Zi);
        n = linspace(0,N-1,N);
        w = 0.5.*(1-cos(2.*pi.*n/(N-1)));
        Zi = Zi.*repmat(w.',1,size(Zi,2));

    case 'humming'
        N = length(Zi);
        n = linspace(0,N-1,N);
        w = 0.54-0.46.*cos(2.*pi.*n/(N-1));
        Zi = Zi.*repmat(w.',1,size(Zi,2));

    case 'blackmann'
       N = length(Zi);
       n = linspace(0,N-1,N); 
       w = 0.35875 - 0.48829.*cos(2.*pi.*n/(N-1))+0.14128.*cos(4.*pi.*n/(N-1))-0.01168.*cos(6.*pi.*n/(N-1));
       Zi = Zi.*repmat(w.',1,size(Zi,2));

    case 'kaiser-bessel'
       N = length(Zi);
       n = linspace(0,N-1,N);
       alpha = 3;
       w = besseli(0,pi.*alpha.*sqrt(1-((2.*n+1-N)./(N-1)).^2))./besseli(0,pi.*alpha);
       Zi = Zi.*repmat(w.',1,size(Zi,2));

end


%% Make FFT

% Former problem of the code when length(Z) = length(T) was uneven number: 
% t_step and thus Fs was computed based on vector T. If length(Zi) was uneven, its last element was removed (see commented code above).
% This led to L = length(fft(Zi)) always being an even number. However, this created a problem when computing vector f: 
% Fs was computed based on T exhibiting an uneven number of elements, while L is the array length AFTER the removal of the last element.
% Therefore, Fs and L did not fit together and the vector f exhibited an error.
% Now, t_step and thus Fs is computed based on the interpolation query points Ti which always has an even number of elements. 
% As a result, Fs and L fit together and f is computed correctly.

    Fs = 1/t_step;                  % Sampling frequency of FFT input signal (can be different from sampling frequency of vector T due to interpolation)

    Y = fft(Zi);
    L = length(Y);
    Y = Y/L;                        
    f = Fs/L*(0:(L/2-1));           % f starts at 0 and goes up to Fs/2 - Fs/L (Fs/2: Nyquist frequency)
   
    S = Y(1:L/2,:);                 % Y(1:L/2+1,:) is the amplitude at the Nyquist frequency
    S(2:end,:) = 2*S(2:end,:);      % Formerly, it was "2:end-1", but end-1 is only correct when f(end) is the Nyquist frequency (not included here)
    
    if strcmpi(output,'polar')
        varargout{1} = abs(S);
        varargout{2} = atan2(-imag(S),real(S));

    elseif  strcmpi(output,'cartesian')
        varargout{1} =  real(S);
        varargout{2} = -imag(S);

    else 
        error('Identifier for output format of the return argument is not known. Possible formats are ''polar'' or ''cartesian''.');
    
    end
    

end