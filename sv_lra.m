function min_R = sv_lra(x,M,display,eps)
%--------------------------------------------------------------------------
% PURPOSE
%  Approximate the 1-dimensional non-uniform discrete Fourier transform 
%  using Low Rank Approximation method. See Ruiz-Antolin and Townsend
%  (2017) for more details.
%
% INPUT: x = [x_1; x_2; ... ;x_N]   non-uniform sampling points in [0,1]
%        M                          number of frequencies w s.t. 
%                                   -M/2 <= w < M/2
%        display                    whether to display the figure or not
%        eps                        precision required
%
% OUTPUT: min_R                     minimum number R required to satisfy 
%                                   the precision eps
%--------------------------------------------------------------------------

%-Parameters---------------------------------------------------------------
% Number of sampling points
N = length(x);

% Set default value
if nargin < 2
    M = N;
end
if nargin < 3
    display = False;
end
if nargin < 4
    eps = 1e-8;
end

%-Generate matrices K and F--------------------------------------------------
% Frequencies
w = (-M/2:M/2-1)';

% Uniform gird
xi = round(x*N)/N;

% Exact NUDFT matrix
F = exp(-2i*pi*w*xi');

% Approximate NUDFT matrix
K = exp(-2i*pi*w*x');

%-Compute the singular values of matrix K./F-------------------------------
% singular values
sv = svd(K./F);

% Minimum number R required
min_R = sum(sv > eps);

%-Display------------------------------------------------------------------
if display
    clf;
    semilogy((1:N),sv,'b-o')
    xlabel('i')
    ylabel('i_{th} singular value')
    title(sprintf((['Singular Values of K./F:  N=%d, eps=%.0e' ...
        '  >>  min R=%d']), N,eps,min_R))
    grid on
    hold on
    yline(eps,'--r')
end
