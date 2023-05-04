function f = gi_method(c,x,M,R,S,tau)
%--------------------------------------------------------------------------
% PURPOSE
%  Approximates the 1-dimensional non-uniform discrete Fourier transform 
%  using Guassian Interpolation method. See Greengard and Lee (2004) for
%  more details.
%
% INPUT: c = [c_1; c_2; ... ;c_N]   1-dimensional input vector(s)
%        x = [x_1; x_2; ... ;x_N]   non-uniform sampling points in [0,1]
%        M                          number of frequencies w s.t. 
%                                   -M/2 <= w < M/2
%        R                          oversampling factor
%        S                          number of grid points to spread data to
%        tau                        spreading factor of Gaussian kernel
%
% OUTPUT: f = [f_1; f_2; ... ;f_M]  Fourier coefficients
%--------------------------------------------------------------------------

%-Parameters---------------------------------------------------------------
% Number of sampling points
N = length(x);

% Number of input vectors
T = size(c,2);

% Set default values
if nargin < 2
    M = N;
end
if nargin < 3
    R = 2;
end
if nargin < 4
    S = 20;
end
if nargin < 5
    tau = (1/M^2)*(pi*(S-1)/2)/(R*(R-0.5));
end

if S == 1
    tau = (1/M^2)*(pi*S/2)/(R*(R-0.5));
end

% Number of points in oversampled grid
M_r = R*M;

% Frequencies
w = (-M/2:M/2-1)';

%-Construct oversampled uniform grid---------------------------------------
% Fine upsampled grid in [0,1]
xi = linspace(0,1,M_r)';

% Step size in upsampled grid
h = 1/M_r;

% Find indices in xi that are the closest to x_j
distances = bsxfun(@(x,y) abs(x-y), x(:), xi');
[~, idx] = min(distances,[],2);

%-Guassian interpolation---------------------------------------------------
f_tau = zeros(M_r,T);
s = fix(S/2);
for j = 1:N
    for l = -s:S-s-1
        f_tau(mod(idx(j)+l,M_r)+1,:) = f_tau(mod(idx(j)+l,M_r)+1,:) + ...
            c(j,:)*exp(-4*pi^2*(x(j)-h*idx(j)-h*l)^2/(4*tau));
    end
end

%-FFT and shift the frequencies -M_r/2 <= k < M_r/2------------------------
F_tau = fft(f_tau,M_r);

%-Correct for spreading in Fourier domain for -M/2 <= k < M/2--------------
f = sqrt(pi/tau) * exp(w.^2*tau) .* F_tau([M_r-M/2+1:M_r,1:M/2],:)/M_r;
