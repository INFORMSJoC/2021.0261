%
% function to numerically compute covariance of EC estimator
% when applying ISDM and G_0 is N(mean0, var0), so
% sum Y ~ N(meansum, varsum).
% 
%
function cov = covisdmnormal(z)

global m im theta beta
global s is ms scale
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

tcdf = 1 - normcdf(z, meansum(is,im), sdsum(is,im) );

fun = @(y) ( y .* normpdf(y, meansum(is,im), sdsum(is,im)) );
    
term1 = ( -tcdf ) .* ( integral(fun, -Inf, z) - meansum(is,im) ...
      .* ( msispar .* normcdf(z, qapp(is,im), sdsum(is,im)) ...
      + (1-msispar) .* normcdf(z, meansum(is,im), sdsum(is,im)) ) );



fun = @(y) ( ( lrisdm(y) - tcdf ) .* ( y .* lrisdm(y) - meansum(is,im) ) ...
    .* ( msispar .* normpdf(y, qapp(is,im), sdsum(is,im)) ...
    + (1-msispar) .* normpdf(y, meansum(is,im), sdsum(is,im))) );

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


