function [D,F,B,e_init,e_ulti] = optimization_nufft(x,M,R,S,init,order, ...
    N_iter,display)
%--------------------------------------------------------------------------
% PURPOSE
%  Use alternating minimization to solve the optimization problem 
%  
%          D, B = argmin || K - D*F*B ||_{Fro}
%                  D, B                                             
%      (K: NUDFT matrix, F: Fourier matrix with permutations)
%  
%  to get the matrix representation of Optimization method.
%
% INPUT: x = [x_1; x_2; ... ;x_N]   non-uniform sampling data in [0,1]
%        M                          number of frequencies w s.t. 
%                                   -M/2 <= w < M/2
%        R                          oversampling factor
%        S                          number of grid points to spread data to
%        init                       initialization method 
%                                   ("gi"(recommended), "id" or "randn")
%        order                      order of alternating minimization 
%                                   ("B-first" or "D-first")
%        N_iter                     number of iterations
%        display                    whether to display the figure or not
%
% OUTPUT: D                         diagonal matrix
%         F                         Fourier matrix (with permutations)
%         B                         sparse matrix
%         e_init                    initial relative error F-norm
%         e_ulti                    ultimate relative error F-norm
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
    S = 12;
end
if nargin < 5
    init = "gi";
end
if nargin < 6
    order = "B_first";
end
if nargin < 7
    N_iter = 5;
end
if nargin < 8
    display = true;
end

w = (-M/2:M/2-1)';
K = exp(-2i*pi*w*x');
M_r = R*M;

%-Record the non-zero entries of sparse matrix B---------------------------
% Construct oversampled uniform grid in [0,1]
xi = linspace(0,1,M_r)';

% Find indices in xi that are the closest to x_j
distances = bsxfun(@(x,y) abs(x-y), x(:), xi');
[~, idx] = min(distances,[],2);

% Record the non-zero entries of matrix B
s = fix(S/2);
idx_list = mod(idx+(-s:S-s-1)-1,M_r) + 1;

%-Optimization-------------------------------------------------------------
% Record the history of error
e = zeros(N_iter+1,1);

% Initialization
[D_gi,F,B_gi,~] = guassian_interpolation(x,M,R,S);

rng(57);

if init == "gi"
    D = D_gi;
    B = B_gi;
elseif init == "eye"
    D = eye(M);
    B = randn(M_r,N).*(B_gi~=0);
    if order == "D-first"
        warning('InputWarning. The setting is equivalent to random initialization.')
    end
elseif init == "randn"
    D = diag(randn(M,1));
    B = randn(M_r,N).*(B_gi~=0);
else
    error('InputError. Variable {init} must be "gi", "eye" or "randn".')
end

e(1) = norm(D*F*B-K,'fro') / sqrt(M*N);

% Iterate using alternating minimization
for n = 1:N_iter
    [D,B] = iter(K,F,D,B,idx_list,order);
    e(n+1) = norm(D*F*B-K,'fro') / sqrt(M*N);
end

%-Record the initial and ultimate relative error norm----------------------
e_init = e(1);
e_ulti = e(n+1);

%-Display------------------------------------------------------------------
if display
    clf;
    figure(1);
    subplot(2,1,1)
    semilogy((0:N_iter),e,'b-o')
    xlabel('number of iterations')
    ylabel('relative error F-norm')
    title(sprintf('Alternating Minimization: iteration 0-%d',N_iter))
    grid on
    subplot(2,1,2)
    semilogy((1:N_iter),e(2:N_iter+1),'b-o')
    xlabel('number of iterations')
    ylabel('relative error F-norm')
    title(sprintf('Alternating Minimization: iteration 1-%d',N_iter))
    grid on

    fprintf(['Initial relative error F-norm:   %.2e\n', ...
             'Ultimate relative error F-norm:  %.2e\n'], ...
             [e(1),e(end)])
end

end


%-Functions----------------------------------------------------------------
% Perform an iteration
function [D,B] = iter(K,F,D,B,idx_list,order)
if order == "B-first"
    B = min_B(K,F,D,idx_list);
    D = min_D(K,F,B);
elseif order == "D-first"
    D = min_D(K,F,B);
    B = min_B(K,F,D,idx_list);
else
    error('InputError. Variable {order} must be "B-first" or "D-first".')
end

end


% Solve matrix B to minimize
function B = min_B(K,F,D,idx_list)

[~,N] = size(K);
[~,M_r] = size(F);
B = zeros(M_r,N);

D_F = D*F;

for i = 1:N
    idxes = idx_list(i,:);
    B(idxes,i) = D_F(:,idxes) \ K(:,i);
end

end


% Solve matrix D to minimize
function D = min_D(K,F,B)

[M,~] = size(K);
D = zeros(M,M);

F_B = F*B;

for i = 1:M
    D(i,i) = K(i,:) / F_B(i,:);
end

end

