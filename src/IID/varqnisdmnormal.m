%
% function to numerically compute variance of quantile estimator
% when applying ISDM and G_0 is N(mean0, var0), so
% sum Y ~ N(meansum, varsum).
% 
%
function var = varqnisdmnormal(z)

global m im theta beta
global s is ms scale
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global msispar dewtq dewtm
global p


tcdf = 1 - normcdf(z, meansum(is,im), sdsum(is,im) );
% temp

% fprintf('\n normcdf(z, meansum(is,im), sdsum(is,im)) = %e', normcdf(z,meansum(is,im), sdsum(is,im)) );

a1 = ( tcdf.^2 ) ...
    .* ( msispar.*normcdf(z,qapp(is,im),sdsum(is,im)) ...
    + (1-msispar).*normcdf(z,meansum(is,im),sdsum(is,im)) );

% fprintf('\n msispar = %e', msispar);
% fprintf('\n theta(is) = %e', theta(is));
% fprintf('\n m(is,%d) = %e', im, m(is,im));
% fprintf('\n cgf0(is) = %e', cgf0(is));
% fprintf('\n normcdf(z, meansum(is,im), sdsum(is,im)) = %e', normcdf(z, meansum(is,im), sdsum(is,im)));
% fprintf('\n normcdf(z, qapp(is,im), sdsum(is,im)) = %e', normcdf(z, qapp(is,im), sdsum(is,im)));



fun = @(y) ( ( ( lrisdm(y) - tcdf ).^2 ) ...
    .* ( msispar .* normpdf(y, qapp(is,im), sdsum(is,im)) ...
    + (1-msispar) .* normpdf(y, meansum(is,im), sdsum(is,im)) ) );
    
a2 = integral(fun, z, Inf);

var = a1 + a2;


end


