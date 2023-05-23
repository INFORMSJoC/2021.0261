%
% function to numerically compute first term in
% exact and approximate variance of importance-sampling 
% estimator of p-quantile of sum of m i.i.d. random variables 
% with marginal CDF G_0 as N(0,1).
% 
%
function a1 = a1normal(z)

global m im theta beta
global s is ms scale
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global p


a1 = ((1.0 - normcdf(z,meansum(is,im), sdsum(is,im))).^2) ...
    .* normcdf(z,qapp(is,im),sdsum(is,im));

end


