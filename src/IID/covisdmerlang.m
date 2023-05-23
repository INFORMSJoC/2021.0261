%
% function to numerically compute covariance of EC estimator
% when applying ISDM and G_0 is Erlang
% 
%
function cov = covisdmerlang(z)

global m im theta beta
global s is ms scale scalenew
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global msispar dewtq dewtm
global p


% fun = @(y) ( y .* exp(m(is,im).*cgf0(is) - theta(is).*y) ...
%     .* normpdf(y, meansum(is,im), sdsum(is,im)) );
%     
% cov = integral(fun, lowlim, Inf) ...
%     - ( 1 - normcdf(lowlim,meansum(is,im),sdsum(is,im)) ) .* meansum(is,im);
% 
% 
% fprintf('\n cov (%d) = %12.5e', m(is,im), cov);
% 


% alternative way of computing covariance

tcdf = 1 - gamcdf( z, ms, scale );

fun = @(y) ( y .* gampdf( y, ms, scale ) );
    
term1 = ( -tcdf ) .* ( integral(fun, -Inf, z) - meansum(is,im) ...
      .* ( msispar .* gamcdf(z, ms, scalenew ) ...
      + (1-msispar) .* gamcdf(z, ms, scale ) ) );



fun = @(y) ( ( lrisdm(y) - tcdf ) .* ( y .* lrisdm(y) - meansum(is,im) ) ...
    .* ( msispar .* gampdf(y, ms, scalenew ) ...
    + (1-msispar) .* gampdf(y, ms, scale ) ) );

term2 = integral(fun, z, Inf);

cov = term1 + term2;

% fprintf('\n cov2(%d) = %12.5e', m(is,im), term1+term2);



% % alternative2 way of computing covariance
% 
% fun = @(y) ( ( (y > lowlim) .* exp( m(is,im).*cgf0(is) - theta(is).*y ) ...
%     - ( 1 - normcdf(lowlim, meansum(is,im), sdsum(is,im)) ) ) ...
%     .* ( y .* exp( m(is,im).*cgf0(is) - theta(is).*y ) - meansum(is,im) ) ...
%     .* normpdf(y, qapp(is,im), sdsum(is,im)) );
% 
% term1 = integral(fun, -Inf, Inf);
% 
% fprintf('\n cov3(%d) = %12.5e', m(is,im), term1);
% 

end


