%
% function to numerically compute variance of quantile estimator
% when applying ISDM and G_0 is Erlang.
% 
%
function var = varqnisdmerlang(z)

global m im theta beta
global s is ms scale scalenew
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global msispar dewtq dewtm
global p


tcdf = 1 - gamcdf(z, ms, scale );
% temp


a1 = ( tcdf.^2 ) ...
    .* ( msispar.*gamcdf( z, ms, scalenew ) ...
    + (1-msispar).*gamcdf( z, ms, scale ) );


fun = @(y) ( ( ( lrisdm(y) - tcdf ).^2 ) ...
    .* ( msispar .* gampdf( y, ms, scalenew ) ...
    + (1-msispar) .* gampdf(y, ms, scale ) ) );
    
a2 = integral(fun, z, Inf);

var = a1 + a2;


end


