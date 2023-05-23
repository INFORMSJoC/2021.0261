%
% function to numerically compute second term in
% exact and approximate variance of importance-sampling 
% estimator of p-quantile of sum of m i.i.d. random variables 
% with marginal CDF G_0 as N(0,1).
% 
%
function a2 = a2normal(z)

global m im theta beta
global s is ms scale
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global p


% fprintf('\n a2normal: normcdf(z, meansum(is,im), sdsum(is,im)) = %e', ...
%     normcdf(z,meansum(is,im), sdsum(is,im)) );


tcdf = 1 - normcdf(z, meansum(is,im), sdsum(is,im));

fun = @(y) ( ( lr(y) - tcdf ).^2 ) ...
    .* normpdf(y, qapp(is,im), sdsum(is,im));
    
a2 = integral(fun, z, Inf);

end
