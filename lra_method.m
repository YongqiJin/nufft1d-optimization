function f = lra_method(c,x,M,R)
%--------------------------------------------------------------------------
% PURPOSE
%  Approximate the 1-dimensional non-uniform discrete Fourier transform 
%  using Low Rank Approximation method. See Ruiz-Antolin and Townsend
%  (2017) for more details.
%
% INPUT: c = [c_1; c_2; ... ;c_N]   1-dimensional input vector(s)
%        x = [x_1; x_2; ... ;x_N]   non-uniform sampling points in [0,1]
%        M                          number of frequencies w s.t. 
%                                   -M/2 <= w < M/2
%        R                          number of items of Taylor expansion
%
% OUTPUT: f = [f_1; f_2; ... ;f_M]  Fourier coefficients
%--------------------------------------------------------------------------

%-Parameters---------------------------------------------------------------
% Number of sampling points
N = length(x);

% Number of input vectors
T = size(c,2);

% Frequencies
w = (-M/2:M/2-1)';

%-Find the closest grid points and record the distances--------------------
h = 1/M;
xi = round(x*M);
idx = mod(xi,M) + 1;
dist = x - xi*h;

%-Compute the Fourier coeffecients-----------------------------------------
f = zeros(M,T);

% Storage the diagnal matrices with vectors
Dx = ones(N,1);
Dk = ones(M,1);

% Compute the Fourier coeffecients using LRA method
for r = 1:R
    v = Dx.*c;
    F = fft(eye(M));
    f = f + Dk.*(F([M/2+1:M,1:M/2],idx)*v);
    Dx = Dx.*(dist.*N);
    Dk = Dk.*(w.*(-2i*pi/N/r));
end
