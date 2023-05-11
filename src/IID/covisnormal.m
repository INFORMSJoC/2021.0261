%
% function to numerically compute covariance of EC estimator
% when applying IS and G_0 is N(mean0, var0), so
% sum Y ~ N(meansum, varsum).
% 
%
function cov = covisnormal(z)

global m im theta beta
global s is ms scale
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
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
% using approach described in "Numerical Issues when
% Computing Asymptotic Variance of IS Quantile Estimator"
% to avoid numerical issues from round-off errors.
% To do this, for z as the true quantile q or its approximation qapp, write 
% gamma(z) 
% = E_{IS} [ ( I( Y > z ) L(Y) - (1-F(z) ) * ( Y L(Y) - mu ) ]
% = b1(z) + b2(z),
% where E_{IS} is the expectation under IS,
% mu is the mean of the iid sum, 
% b1(z) 
% = E_{IS} [ ( I( Y > z ) L(Y) - (1-F(z) ) * ( Y L(Y) - mu ) ; Y <= z ],
% and
% b2(z) 
% = E_{IS} [ ( I( Y > z ) L(Y) - (1-F(z) ) * ( Y L(Y) - mu ) ; Y > z ].
% After some simplifications, we can write
% b1(z) 
% = -(1-F(z)) * int_{-\infty}^{z} y dF(y) + (1-F(z)) * mu * F_{IS}(z)
% and
% b2(z)
% = int_{z}^{\infty} (L(y) - (1-F(z))*(y*L(y) - mu) dF_{IS}(y),
% where the integrals in both b1(z) and b2(z) have nonnegative
% integrands.


tcdf = 1 - normcdf( z, meansum(is,im), sdsum(is,im) );

fun = @(y) ( y .* normpdf(y, meansum(is,im), sdsum(is,im)) );
    
% fun = @(y) ( y .* exp( m(is,im).*cgf0(is) - theta(is).*y ) - meansum(is,im) )...
%     .* normpdf(y, qapp(is,im), sdsum(is,im));

term1 = ( -tcdf ) .* ( integral(fun, -Inf, z) ...
      - meansum(is,im) .* normcdf(z, qapp(is,im), sdsum(is,im)) );


fun = @(y) ( ( lr(y) - tcdf ) .* ( y .* lr(y) - meansum(is,im) ) ...
    .* normpdf(y, qapp(is,im), sdsum(is,im)) );

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


