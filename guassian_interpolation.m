function [D,F,B,e] = guassian_interpolation(x,M,R,S,tau)
%--------------------------------------------------------------------------
% PURPOSE
%  Generate the matrix representation of Guassian Interpolation method.
%
% INPUT: x = [x_1; x_2; ... ;x_N]   non-uniform positions in [0,1]
%        M                          number of frequencies w s.t. 
%                                   -M/2 <= w < M/2
%        R                          oversampling factor
%        S                          number of grid points to spread data to
%        tau                        spreading factor of Gaussian kernel
%
% OUTPUT: D                         diagonal matrix
%         F                         Fourier matrix (with permutations)
%         B                         sparse interpolation matrix
%         e                         relative error F-norm
%--------------------------------------------------------------------------

%-Parameters---------------------------------------------------------------
% Number of data points
N = length(x);

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
% Construct oversampled uniform grid in [0,1]
xi = linspace(0,1,M_r)';

% Step size in upsampled grid
h = 1/M_r;

% Find indices in xi that are the closest to x_j
distances = bsxfun(@(x,y) abs(x-y), x(:), xi');
[~, idx] = min(distances,[],2);

%-Guassian Interpolation---------------------------------------------------
B = zeros(M_r,N);
s = fix(S/2);
for j = 1:N
    for l = -s:S-s-1
        B(mod(idx(j)+l,M_r)+1,j) = exp(-4*pi^2*(x(j)-h*idx(j)-h*l)^2/(4*tau));
    end
end

%-FFT and shift the frequencies -M_r/2 <= k < M_r/2------------------------
F = fft(eye(M_r));
F = F([M_r-M/2+1:M_r,1:M/2],:);

%-Correct for spreading in Fourier domain for -M/2 <= k < M/2--------------
D = diag( sqrt(pi/tau) .* exp(w.^2*tau) / M_r);

%-Compute exact NUDFT matrix-----------------------------------------------
K = exp(-2i*pi*w*x');

%-Compute relative error F-norm--------------------------------------------
e = norm(D*F*B - K, 'fro') / sqrt(M*N);
