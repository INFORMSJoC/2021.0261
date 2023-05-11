%
% function to numerically compute covariance of EC estimator
% when applying SRS and G_0 is N(mean0, var0), so
% sum Y ~ N(meansum, varsum).
% 
%
function cov = covsrsnormal(uplim)

global m im theta beta
global s is ms scale
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global p


fun = @(y) ( y .* normpdf(y, meansum(is,im), sdsum(is,im)) );
    
cov = integral(fun, uplim, Inf ) - ( 1 - p(is,im) ) .* meansum(is,im);

end


