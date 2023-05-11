%
% function to numerically compute first term in
% exact and approximate variance of importance-sampling 
% estimator of p-quantile of sum of m i.i.d. random variables 
% with marginal CDF G_0 as Erlang with.
% 
%
function a1 = a1erlang(z)

global m im theta beta
global meansum sdsum
global s is ms scale scalenew
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global p


a1 = ((1.0 - gamcdf(z,ms,scale)).^2) .* gamcdf(z,ms,scalenew);
% fprintf('\n ms = %d, gamcdf(z,ms,scale) = %e, gamcdf(z,ms,scalenew) = %e', ...
%     ms, gamcdf(z,ms,scale), gamcdf(z,ms,scalenew));

end