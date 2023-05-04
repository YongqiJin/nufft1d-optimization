function e = mean_relative_error_norm(f_hat,f)
%--------------------------------------------------------------------------
% PURPOSE
%  Compute the mean relative error L2-norm between x_hat and x.
%
% INPUT: x_hat                      calculated vector(s)
%        x                          exact vector(s)
%
% OUTPUT: e                         the mean relative L2 norm
%--------------------------------------------------------------------------

%-Parameters---------------------------------------------------------------
% Number of vectors
T = size(f,2);

%-Compute the mean relative error L2-norm----------------------------------
e = zeros(1,T);

for i = 1:T
    e(i) = norm(f_hat(:,i)-f(:,i),2) / norm(f(:,i),2);
end

e = mean(e);
