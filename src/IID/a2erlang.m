%
% function to numerically compute second term in
% exact and approximate variance of importance-sampling 
% estimator of p-quantile of sum of m i.i.d. random variables 
% with marginal CDF G_0 as Erlang.
% 
%
function a2 = a2erlang(z)

global m im theta beta
global meansum sdsum
global s is ms scale scalenew
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global p


fun = @(y) ( ( exp(m(is,im).*cgf0(is) - theta(is).*y) ...
    - 1 + gamcdf(z,ms,scale) ).^2 ) ...
    .* gampdf(y, ms, scalenew);

a2 = integral(fun, z, Inf);


end
