function f = opt_method(c,D,F,B)
%--------------------------------------------------------------------------
% PURPOSE
%  Approximate the 1-dimensional non-uniform discrete Fourier transform 
%  using Optimization method. 
%
% INPUT: c = [c_1; c_2; ... ;c_N]   1-dimensional input vector(s)
%        D                          diagonal matrix
%        F                          Fourier matrix (with permutations)
%        B                          sparse matrix
%
% OUTPUT: f = [f_1; f_2; ... ;f_M]  Fourier coefficients
%--------------------------------------------------------------------------

f = D*((F*(B*c)));
