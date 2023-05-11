%
% function to numerically compute variance of mean estimator
% when applying ISDM and G_0 is N(mean0, var0), so
% sum Y ~ N(meansum, varsum).
% 
%
function var = varmisdmnormal

global m im theta beta
global s is ms scale
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global msispar dewtq dewtm
global p


% fun = @(y) ( (y.^2) .* lrisdm(y)  ...
%     .* normpdf(y, meansum(is,im), sdsum(is,im)) );
%     
% var = integral(fun, -Inf, Inf) - meansum(is,im).^2;
% 
% fprintf('\n varmisdsm (%d) = %e', im, var);
% 
% 


fun = @(y) ( ( ( y .* lrisdm(y) - meansum(is,im) ).^2 ) ...
    .* ( msispar .* normpdf(y, qapp(is,im), sdsum(is,im)) ...
    + (1-msispar) .* normpdf(y, meansum(is,im), sdsum(is,im)) ) );
    
var = integral(fun, -Inf, Inf);

% fprintf('\n varmisdsm2(%d) = %e', im, var);

end


