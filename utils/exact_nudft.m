function f = exact_nudft(c,x,M)
%--------------------------------------------------------------------------
% PURPOSE
%  Compute the discrete Fourier transform of the uniform or non-uniform
%  sampling data c exactly.
%
% INPUT: c = [c_1; c_2; ... ;c_N]   1-dimensional input vector(s)
%        x = [x_1; x_2; ... ;x_N]   uniform or non-uniform positions 
%                                   in [0,1]
%        M                          number of frequencies w s.t.
%                                   -M/2 <= w < M/2
%
% OUTPUT: f = [f_1; f_2; ... ;f_M]  Fourier coefficients
%--------------------------------------------------------------------------

%-Parameters---------------------------------------------------------------
% Number of data points
N = length(c);

% Set default values
if nargin < 2
    x = (0:N-1)'/N;
    x = 2*pi*x;
end
if nargin < 3
    M = N;
end

% Frequencies
w = (-M/2:M/2-1)';

%-Compute the Fourier coefficients-----------------------------------------
% exact NUDFT matrix
K = exp(-2i*pi*w*x');

% Fourier coefficients
f = K*c;
