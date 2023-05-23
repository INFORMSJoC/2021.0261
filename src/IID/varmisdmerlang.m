%
% function to numerically compute variance of mean estimator
% when applying ISDM and G_0 is Erlang
% 
%
function var = varmisdmerlang

global m im theta beta
global s is ms scale scalenew
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global msispar dewtq dewtm
global p



fun = @(y) ( ( ( y .* lrisdm(y) - meansum(is,im) ).^2 ) ...
    .* ( msispar .* gampdf(y, ms, scalenew ) ...
    + (1-msispar) .* gampdf(y, ms, scale ) ) );
    
var = integral(fun, -Inf, Inf);
% var = integral(fun, -q(is,im), q(is,im));

% fprintf('\n varmisdsm2(%d) = %e', im, var);

end


