%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function main:
%
% written by Marvin Nakayama
%
% Code for generating numerical results for the i.i.d. sum model
% in the paper, 
% Li, Kaplan, and Nakayama, "Monte Carlo Methods for Economic Capital",
% to appear in INFORMS Journal on Computing.
%
% Overview:
% Let Y ~ F be a sum of m i.i.d. random variables, each with
% marginal CDF G_0. The code numerically computes (using exact 
% analytical formulas or quadrature) the following quantities:
% mean of Y,
% p-quantile of Y,
% ec (economic capital) of Y.
% The quantile level is defined as p = 1 - exp(-beta*m),
% where beta > 0 is a constant specified in the code, as in Glynn (1996).
% Also, the code numerically computes the asymptotic variance of the
% estimators of each of the above for various Monte Carlo (MC) methods:
% simple random sampling (srs),
% importance sampling (is) with exponential twist based on Glynn (1996),
% measure-specific importance sampling (msis),
% importance sampling with defensive mixture (isdm),
% double estimator (de) with fixed weights or optimal weights (deo).
% From these, the code computes the relative error (RE) of each
% estimator. The output is not from simulation experiments but
% rather through numerical computations.
%
% The code also computes various approximations for the p-quantile, 
% the ec, and the asymptotic variances and relative errors
% of their estimators for each of the above MC methods.
% One approximation replaces the true p-quantile (q) with
% the quantile approximation (often denoted qapp) given by
% qapp = m * Q_0(theta),
% where m is the number of i.i.d. summands,
% Q_0 is the cumulative generating function of a single summand,
% and theta is the root of the equation
% -theta * Q_0'(theta) + Q_0(theta) = -beta,
% as proposed by Glynn (1996) for estimating the
% p-quantile using exponential twisting with parameter theta.
% Another approximation (fapp) uses a saddlepoint approximation 
% to the true density f of the i.i.d. sum Y.
% The two approximations are sometimes combined as fappqapp,
% which is the approximate density evaluated at the approximate quantile.
%
% The code works for two choices of the marginal distribution G_0
% of each summand:
% N(mean0,var0): normal with mean mean0 = 1, variance var0 = 1,
% and 
% Gamma(s, 1): shape parameter = s, scale parameter = 1,
% where scale parameter is the mean of each stage of the Erlang.
% The quantile level p = 1 - exp(-beta.*m(is,im)),
% where beta > 0 is a specified constant.
% 
% The code works for values of beta only below 
% some threshold depending on parameters of the model, 
% including the marginal CDF G_0, its parameters, and the number m
% of i.i.d. summands. When beta or m is too large, the program 
% has numerical issues leading to NaN outputs arising 
% from overflow/underflow because of the likelihood ratio 
% growing too large/small.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% m: array of numbers of iid summands, in powers of 2
global m 

% im: index into array m
global im 

% immax: size of array m
global immax 

% theta: IS twisting parameter, computed as positive root of
% -theta * Q_0'(theta) + Q_0(theta) = -beta, for fixed beta > 0
global theta 

% beta: used in computing quantile level p = 1 - exp(-beta*m)
global beta

% meansum: exact mean of iid sum
global meansum 

% sdsum: standard deviation of iid sum
global sdsum

% varqsrs: CLT asymptotic variance of SRS quantile estimator
global varqsrs

% s: matrix of numbers of stages in Erlang
global s 

% is: index into matrix s
global is 

% ms: product of m and s (used only for Erlang)
global ms 

% scale: original scale parameter of Erlang (i.e., mean of single stage)
global scale 

% scalenew: scale parameter of Erlang under IS with exponential twist
global scalenew

% mgf0: moment generating function of marginal distn G_0 of single summand
global mgf0 

% cgf0: cumulant generating function of G_0
global cgf0 

% cgf0d: derivative of cumulant generating function of G_0
global cgf0d 

% cgf0d2: second derivative of cumulant generating function of G_0
global cgf0d2 

% q: true value of p-quantile
global q 

% qapp: approximate value of p-quantile, computed as qapp = m * cgf0d(theta)
global qapp

% fq: density of sum evaluated at true quantile
global fq 

% fappqapp: approximate density of sum evaluated at approximate quantile
global fappqpp

% msispar: sampling allocation parameter in [0,1] for measure-specific IS
global msispar 

% dewtq: double estimator weight combining IS and SRS quantile estimators
global dewtq 

% dewtm: double estimator weight combining IS and SRS mean estimators
global dewtm

% p: quantile level
global p

% idist: indicator for G_0, either = inormal, or = ierlang
global idist 
global inormal 
global ierlang

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameter for quantile level p = p(is,im).
%
% beta = 1;
beta = 1.1;
% beta = 2.2;
% beta = 3.3;
% beta = 6.6029;
% beta = 7.0;
% beta = 7.1;

% immax: number of different values of m (powers of 2) to compute,
% where m is the number of iid summands.
% 
immax = 6;

% msis sampling allocation parameter, where for a given
% total sample size of n, msispar * n is used for IS,
% and (1-msispar)*n is used for SRS. 
% Also used for isdm.
%
msispar = 0.5;

% de weights for combining quantile estimators and mean estimators.
% The the DE quantile estimator is a convex combination of
% the IS quantile estimator (with weight dewtq) and 
% the SRS quantile estimator (with weight (1-dewtq)).
% The the DE mean estimator is a convex combination of
% the IS mean estimator (with weight dewtm) and 
% the SRS mean estimator (with weight (1-dewtm)).
% 
dewtq = 0.5;
dewtm = 0.5;

% indicator for marginal distribution of each summand
inormal = 1;
ierlang = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idist = inormal;

% is: used mainly for Erlang, so set is = 1 for normal
is = 1;

% specify mean and variance of marginal distn of each summand
% mean0 = 0;
mean0 = 1;
var0 = 1;
sd0 = sqrt(var0);

% twisting parameter for IS using Glynn (1996) for IS quantile estimation
theta(is) = sqrt(2*beta)/sd0;

% MGF(theta) for G_0
mgf0(is) = exp( ( mean0.*theta(is) ) + ( var0.*(theta(is).^2)/2 ) );

% MGF(-theta) for G_0, first derivative, and second derivative
% (the trailing "n" means negative argument)
%
mgf0n(is) = exp( ( -mean0.*theta(is) ) + ( var0.*(theta(is).^2)/2 ) );
mgf0dn(is) = ( mean0 - var0.*theta(is) ) .* mgf0n(is);
mgf0d2n(is) = var0 .* mgf0n(is) ...
    + ( ( mean0 - var0.*theta(is) ).^2 ) .* mgf0n(is);

% CGF(theta) for G_0 and first two derivatives
cgf0(is) = ( mean0.*theta(is) ) + ( var0.*(theta(is).^2)/2 );
cgf0d(is) = mean0 + var0.*theta(is);
cgf0d2(is) = var0;

% CGF(-theta) for G_0 and derivatives
% (the trailing "n" means negative argument)
%
cgf0n(is) = ( -mean0.*theta(is) ) + ( var0.*(theta(is).^2)/2 );
cgf0dn(is) = mean0 - var0.*theta(is);
cgf0d2n(is) = var0;

% fprintf('alpha = %9.3e', mgf0(is).*mgf0n(is));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normal: RE computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop through different number of i.i.d. summands
for im = 1:immax
    
    % compute number m(is,im) of summands as power of 2, 
    % and also compute quantile level p as in Glynn (1996)
    %
    m(is,im) = 2.^(im-1);
    p(is,im) = 1.0 - exp(-beta.*m(is,im));

    % compute mean and variance of i.i.d. sum
    %
    meansum(is,im) = m(is,im).*mean0;
    varsum(is,im) = m(is,im).*var0;
    sdsum(is,im) = sqrt(varsum(is,im));
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% q: srs (simple random sampling) exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute exact quantile and exact density at exact quantile
    q(is,im) = norminv(p(is,im), meansum(is,im), sdsum(is,im));
    fq(is,im) = normpdf(q(is,im), meansum(is,im), sdsum(is,im));

    % compute exact asymptotic variance of SRS quantile estimator
    varqsrs(is,im) = p(is,im) .* ( 1 - p(is,im) ) / ( fq(is,im).^2 );
    sdqsrs(is,im) = sqrt( varqsrs(is,im) );
    reqsrs(is,im) = sdqsrs(is,im) / q(is,im);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% q: srs approx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute quantile approximation based on Glynn (1996)
    qapp(is,im) = m(is,im).*cgf0d(is);
    
    % compute exact density at approximate quantile
    fqapp(is,im) = normpdf(qapp(is,im), meansum(is,im), sdsum(is,im));
    
    % compute saddlepoint approx to density at approximate quantile
    fappqapp(is,im) = ( 1 / sqrt( 2.*pi.*m(is,im).*cgf0d2(is) ) ) ...
        .* exp(m(is,im).*(cgf0(is) - theta(is).*cgf0d(is)));
    
%     fprintf('\n m(is,im) = %d, cgf0 = %e, cgf0d = %e, cgf0d2  = %e', ...
%         m(is,im), cgf0(is), cgf0d(is), cgf0d2(is) );
%     fprintf('\n m(is,im) = %d, fappqapp  = %e', m(is,im), fappqapp(is,im));
%     temp = ( 1 / sqrt( 2.*pi.*m(is,im).*cgf0d2(is) ) ) ...
%         .* exp(- beta.*m(is,im));
%     fprintf('\n m(is,im) = %d, fappqapp2 = %e', m(is,im), temp);


    % approx asymptotic variance of SRS estimator of quantile
    varqsrsapp(is,im) = p(is,im).*(1-p(is,im))/( fappqapp(is,im).^2 );
    sdqsrsapp(is,im) = sqrt( varqsrsapp(is,im) );
    reqsrsapp(is,im) = sdqsrsapp(is,im)/qapp(is,im);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% q: is (importance sampling) exact
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerically computing the asymptotic variance and covariance
% when applying IS can lead to numerical issues, especially when 
% m or beta is large. Because of this, we have to use formulas
% for the variance and covariance that avoids subtractions,
% which can lead to round-off errors. Instead, the code below
% computes these quantities by numerically integrating only
% nonnegative integrands. For example, the numerator of the
% asymptotic variance for the IS quantile estimator is
% num = E_{IS} [ L^2 I(Y > q) ] - (1 - p)^2, 
% where E_{IS} denotes expectation under IS, 
% L = L(Y) is the likelihood ratio, and Y is the iid sum.
% The numerator num is also a variance, but taking the difference 
% as above can lead to negative values because of round-off errors.
%
% Instead, write the numerator in an algebraically equivalent form as
% num = E_{IS} [ ( L(Y) I(Y > q) - (1-p) )^2 ],
% which is an integral whose integrand is always nonnegative,
% so it should lead to num evaluating to something nonnegative.
% To handle the numerator in a similar manner for both the 
% true quantile q and the approximate quantile qapp, instead write
% num(z) = E_{IS} [ ( L(Y) I(Y > z) - (1-p) )^2 ],
% where we set z = q or z = qapp.
% We then write the integral num(z) as the sum of two integrals
% over different ranges for Y, i.e.,
% num(z) = a1(z) + a2(z), where
% a1(z) = E_{IS} [ ( L(Y) I(Y > z) - (1-p) )^2 ; Y <= z ]
% and
% a2(z) = E_{IS} [ ( L(Y) I(Y > z) - (1-p) )^2 ; Y > z ].
% After some simplification, we can write
% a1(z) = ( 1 - F(z) )^2 * F_{IS}(z),
% a2(z) = int_{y=z}^{\infty} ( L(y) - (1-F(z)) )^2 dF_{IS}(y),
% where F_{IS} is the CDF of Y under IS.
% More details appear in our document, 
% "Numerical Computation of Asymptotic Variance of 
% IS Quantile Estimator and Covariance".
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    z = q(is,im);
    
    a1z = a1normal(z);
    a2z = a2normal(z);
    
    varqnum(is,im) = a1z + a2z;
    sdqnum(is,im) = sqrt(varqnum(is,im));
    varqis(is,im) = varqnum(is,im)/( fq(is,im).^2 );
    sdqis(is,im) = sqrt( varqis(is,im) );
    reqis(is,im) = sdqis(is,im)/abs(q(is,im));

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Normal: RE   
% q: is approximations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    z = qapp(is,im);
    
    a1z = a1normal(z);    
    a2z = a2normal(z);
    
    varqnumapp(is,im) = a1z + a2z;
    sdqnumapp(is,im) = sqrt(varqnumapp(is,im));
    varqisapp(is,im) = varqnumapp(is,im)/( fappqapp(is,im).^2 );
    sdqisapp(is,im) = sqrt( varqisapp(is,im) );
    reqisapp(is,im) = sdqisapp(is,im)/abs(qapp(is,im));

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: old formula for varqnum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The next few lines use numerically unstable formulas
% for computing the asymptotic variance of the quantile estimator.
% These are not used, but instead the quantities
% are computed using the more numerically stable
% code above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    m2 = exp(m(is,im).*(cgf0(is)+(theta(is).^2)/2)) .* (1 - normcdf(z,-m(is,im).*theta(is),sdsum(is,im)));
    m1 = 1 - normcdf(z, meansum(is,im), sdsum(is,im));
    
    varnum2 = m2 - (m1.^2);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% mean: srs (simple random sampling) exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    varmsrs(is,im) = varsum(is,im);
    
    sdmsrs(is,im) = sqrt( varmsrs(is,im) );
    remsrs(is,im) = sdmsrs(is,im)/abs(meansum(is,im));
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% mean: is exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (m(is,im) == 1)    
        varmis(is,im) = mgf0(is).*mgf0d2n(is) - mean0.^2;
    else
        alpha = mgf0(is).*mgf0n(is);
        varmis(is,im) = m(is,im) .* (alpha.^m(is,im)) ...
            .* ( (m(is,im) .* (cgf0dn(is).^2)) + cgf0d2n(is) ) ...
            - ( meansum(is,im).^2 );
    end
        
    sdmis(is,im) = sqrt(varmis(is,im));
    remis(is,im) = sdmis(is,im)/abs(meansum(is,im));

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% ec: srs (simple random sampling) exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ec(is,im) = q(is,im) - meansum(is,im);
    
    % covecsrs = E[ I(Y > q) Y ] - (1-p) \mu
    % i.e., not divided by f(q)
    covecsrs(is,im) = covsrsnormal(q(is,im));
    
%     fprintf('\n m = %d, covecsrs = %e', ...
%         m(is,im), covecsrs(is,im) );
    
    varecsrs(is,im) = varqsrs(is,im) + varmsrs(is,im) ...
        - 2.*covecsrs(is,im)/fq(is,im);
    sdecsrs(is,im) = sqrt( varecsrs(is,im) );
    reecsrs(is,im) = sdecsrs(is,im)/abs(ec(is,im));
    

    corrsrs(is,im) = (covecsrs(is,im)/fq(is,im)) ./ ...
    sqrt(varqsrs(is,im).*varmsrs(is,im));

%     fprintf('\n normal, m = %d, corrsrs = %e, corrsrs^2 = %e', ...
%         m(is,im), corrsrs(is,im),  corrsrs(is,im).^2);
%     fprintf('\n');
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% ec: srs approx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ecapp(is,im) = qapp(is,im) - meansum(is,im);
    covecsrsapp(is,im) = covsrsnormal(q(is,im));
    
    varecsrsapp(is,im) = varqsrsapp(is,im) + varmsrs(is,im) ...
        - 2.*covecsrsapp(is,im)/fappqapp(is,im);
    sdecsrsapp(is,im) = sqrt( varecsrsapp(is,im) );
    reecsrsapp(is,im) = sdecsrsapp(is,im)/abs(ecapp(is,im));
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% ec: is (importance sampling) exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % covecis = E[ I(Y > q) Y L ] - (1-p) \mu
    % i.e., not divided by f(q)
    covecis(is,im) = covisnormal(q(is,im));
    
    varecis(is,im) = varqis(is,im) + varmis(is,im) ...
        - 2.*covecis(is,im)/fq(is,im);

%     fprintf('\n varecis (%d) = %10.4e', m(is,im), varecis(is,im));
%     fprintf('\n varecis2(%d) = %10.4e', m(is,im), varecisnormal);
    
    sdecis(is,im) = sqrt(varecis(is,im));
    reecis(is,im) = sdecis(is,im)/abs(ec(is,im));

    
    corris(is,im) = (covecis(is,im)./fq(is,im))./ ...
    sqrt(varqis(is,im).*varmis(is,im));

%     fprintf('\n normal, m = %d, corris = %e, corris^2 = %e', ...
%         m(is,im), corris(is,im),  corris(is,im).^2);
%     fprintf('\n');
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% ec: is approx 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    covecisapp(is,im) = covisnormal(qapp(is,im));
    
    varecisapp(is,im) = varqisapp(is,im) + varmis(is,im) ...
        - 2.*covecisapp(is,im)/fappqapp(is,im);
    
    sdecisapp(is,im) = sqrt(varecisapp(is,im));
    reecisapp(is,im) = sdecisapp(is,im)/abs(ecapp(is,im));

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% ec: msis exact (measure-specific IS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    varecmsis(is,im) = varqis(is,im)./msispar ...
        + varmsrs(is,im)./(1-msispar);
    sdecmsis(is,im) = sqrt( varecmsis(is,im) );
    reecmsis(is,im) = sdecmsis(is,im)/abs(ec(is,im));
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% ec: msis approx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    varecmsisapp(is,im) = varqisapp(is,im)/msispar ...
        + varmsrs(is,im)/(1-msispar);
    sdecmsisapp(is,im) = sqrt( varecmsisapp(is,im) );
    reecmsisapp(is,im) = sdecmsisapp(is,im)/abs(ecapp(is,im));
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% ec: isdm exact (IS with defensive mixture)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    z = q(is,im);
    den = fq(is,im);
    
    varqisdm(is,im) = varqnisdmnormal(z) ./ (den.^2);
    varmisdm(is,im) = varmisdmnormal;
    covisdm(is,im) = covisdmnormal(z);

    varecisdm(is,im) = varqisdm(is,im) + varmisdm(is,im) ...
        - 2.*covisdm(is,im) ./ den;
    sdecisdm(is,im) = sqrt( varecisdm(is,im) );
    reecisdm(is,im) = sdecisdm(is,im)/ec(is,im);
       
    sdqisdm(is,im) = sqrt( varqisdm(is,im) );
    reqisdm(is,im) = sdqisdm(is,im)/abs(q(is,im));
    
    sdmisdm(is,im) = sqrt( varmisdm(is,im) );
    remisdm(is,im) = sdmisdm(is,im)/abs(meansum(is,im));
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% ec: isdm approx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    z = qapp(is,im);
    den = fappqapp(is,im);
    
    varqisdmapp(is,im) = varqnisdmnormal(z) ./ (den.^2);
    varmisdmapp(is,im) = varmisdmnormal;
    covisdmapp(is,im) = covisdmnormal(z);
    
    varecisdmapp(is,im) = varqisdmapp(is,im) + varmisdmapp(is,im) ...
        - 2.*covisdmapp(is,im) ./ den;
    sdecisdmapp(is,im) = sqrt( varecisdmapp(is,im) );
    reecisdmapp(is,im) = sdecisdmapp(is,im)/ecapp(is,im);
       
    sdqisdmapp(is,im) = sqrt( varqisdmapp(is,im) );
    reqisdmapp(is,im) = sdqisdmapp(is,im)/abs(qapp(is,im));
        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% ec: de exact (double estimator with fixed weights)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    z = q(is,im);
    den = fq(is,im);
    
    temp1 = (1 ./ msispar ) ...
        .* ( (dewtq.^2) .* varqis(is,im) ...
        + (dewtm.^2) .* varmis(is,im) ...
        - 2 .* dewtq .* dewtm .* covecis(is,im) ./ den );
    
    temp2 = (1 ./ (1 - msispar) ) ...
        .* ( ((1-dewtq).^2) .* varqsrs(is,im) ...
        + ((1-dewtm).^2) .* varmsrs(is,im) ...
        - 2 .* (1-dewtq) .* (1-dewtm) .* covecsrs(is,im) ./ den );
    
    varqde(is,im) = (1 ./ msispar ) .* (dewtq.^2) .* varqis(is,im) ...
        + (1 ./ (1 - msispar) ) .* ((1-dewtq).^2) .* varqsrs(is,im);
    
    varmde(is,im) = (1 ./ msispar ) .* (dewtm.^2) .* varmis(is,im) ...
        + (1 ./ (1 - msispar) ) .* ((1-dewtm).^2) .* varmsrs(is,im);
    
    covde(is,im) = (1 ./ msispar ).* dewtq .* dewtm .* covecis(is,im) ...
        + (1 ./ (1 - msispar) ) .* (1-dewtq) .* (1-dewtm) .* covecsrs(is,im);
    
    varecde(is,im) = varqde(is,im) + varmde(is,im) ...
        - 2 .* covde(is,im) / den;
    sdecde(is,im) = sqrt( varecde(is,im) );
    reecde(is,im) = sdecde(is,im)/abs(ec(is,im));
       
    sdqde(is,im) = sqrt( varqde(is,im) );
    reqde(is,im) = sdqde(is,im)/abs(q(is,im));
    
    sdmde(is,im) = sqrt( varmde(is,im) );
    remde(is,im) = sdmde(is,im)/abs(meansum(is,im));
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% ec: de approx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    z = qapp(is,im);
    den = fappqapp(is,im);
        
    temp1 = (1 ./ msispar ) ...
        .* ( (dewtq.^2) .* varqisapp(is,im) ...
        + (dewtm.^2) .* varmis(is,im) ...
        - 2 .* dewtq .* dewtm .* covecisapp(is,im) ./ den );
    
    temp2 = (1 ./ (1 - msispar) ) ...
        .* ( ((1-dewtq).^2) .* varqsrsapp(is,im) ...
        + ((1-dewtm).^2) .* varmsrs(is,im) ...
        - 2 .* (1-dewtq) .* (1-dewtm) .* covecsrsapp(is,im) ./ den );
    
    varqdeapp(is,im) = (1 ./ msispar ) .* (dewtq.^2) .* varqisapp(is,im) ...
        + (1 ./ (1 - msispar) ) .* ((1-dewtq).^2) .* varqsrsapp(is,im);
    
    varmdeapp(is,im) = (1 ./ msispar ) .* (dewtm.^2) .* varmis(is,im) ...
        + (1 ./ (1 - msispar) ) .* ((1-dewtm).^2) .* varmsrs(is,im);
    
    covdeapp(is,im) = (1 ./ msispar ).* dewtq .* dewtm .* covecisapp(is,im) ...
        + (1 ./ (1 - msispar) ) .* (1-dewtq) .* (1-dewtm) .* covecsrsapp(is,im);
    
    varecdeapp(is,im) = varqdeapp(is,im) + varmdeapp(is,im) ...
        - 2 .* covdeapp(is,im) ./ den;
    sdecdeapp(is,im) = sqrt( varecdeapp(is,im) );
    reecdeapp(is,im) = sdecdeapp(is,im)/abs(ecapp(is,im));
      
    sdqdeapp(is,im) = sqrt( varqdeapp(is,im) );
    reqdeapp(is,im) = sdqdeapp(is,im)/abs(qapp(is,im));
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal: RE   
% ec: de optimal exact (double estimator with optimal weights)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    z = q(is,im);
    den = fq(is,im);
    
    vqis = varqis(is,im) ./ (msispar);
    vqsrs = varqsrs(is,im) ./ ((1 - msispar));
    vmis = varmis(is,im) ./ msispar;
    vmsrs = varmsrs(is,im) ./ (1 - msispar);
    cis = covecis(is,im) ./ (msispar .* den);
    csrs = covecsrs(is,im) ./ ((1 - msispar) .* den);
    
    dea0 = vqsrs.*vmis - cis.*cis - 2.*cis.*csrs - csrs.*csrs ...
        + vqis.*vmis + vqis.*vmsrs + vqsrs.*vmsrs;
    
    dea1 = vqsrs.*vmis + vqsrs.*vmsrs - vmis.*csrs + vmsrs.*cis ...
        - cis.*csrs - csrs.*csrs;
    
    dea2 = vqis.*vmsrs + vqsrs.*vmsrs - vqis.*csrs + vqsrs.*cis ...
        - cis.*csrs - csrs.*csrs;
    
    dewtqo = dea1 ./ dea0;
    dewtmo = dea2 ./ dea0;
    
%     temp1 = (1 ./ msispar ) ...
%         .* ( (dewtqo.^2) .* varqis(is,im) ...
%         + (dewtmo.^2) .* varmis(is,im) ...
%         - 2 .* dewtqo .* dewtmo .* covecis(is,im) ./ den );
%     
%     temp2 = (1 ./ (1 - msispar) ) ...
%         .* ( ((1-dewtqo).^2) .* varqsrs(is,im) ...
%         + ((1-dewtmo).^2) .* varmsrs(is,im) ...
%         - 2 .* (1-dewtqo) .* (1-dewtmo) .* covecsrs(is,im) ./ den );
    
    varqdeo(is,im) = (1 ./ msispar ) .* (dewtqo.^2) .* varqis(is,im) ...
        + (1 ./ (1 - msispar) ) .* ((1-dewtqo).^2) .* varqsrs(is,im);
    
    varmdeo(is,im) = (1 ./ msispar ) .* (dewtmo.^2) .* varmis(is,im) ...
        + (1 ./ (1 - msispar) ) .* ((1-dewtmo).^2) .* varmsrs(is,im);
    
    covdeo(is,im) = (1 ./ msispar ).* dewtqo .* dewtmo .* covecis(is,im) ...
        + (1 ./ (1 - msispar) ) .* (1-dewtqo) .* (1-dewtmo) .* covecsrs(is,im);
    
    varecdeo(is,im) = varqdeo(is,im) + varmdeo(is,im) ...
        - 2 .* covdeo(is,im) / den;
    sdecdeo(is,im) = sqrt( varecdeo(is,im) );
    reecdeo(is,im) = sdecdeo(is,im)/abs(ec(is,im));
       
    sdqdeo(is,im) = sqrt( varqdeo(is,im) );
    reqdeo(is,im) = sdqdeo(is,im)/abs(q(is,im));
    
    sdmdeo(is,im) = sqrt( varmdeo(is,im) );
    remdeo(is,im) = sdmdeo(is,im)/abs(meansum(is,im));
    
    

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal: RE output: q
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n q: Method1 = srs, Method2 = is');
fprintf('\n G_0 = normal(%d,%d): theta(is) = %10.4e, beta = %10.4e', ...
    mean0, var0, theta(is), beta);

outrehead;

for im = 1:immax
    
    fprintf('\n %2d & %8.2e', m(is,im), 1-p(is,im));
    fprintf(' & %8.2e & %8.2e', q(is,im), qapp(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', reqsrs(is,im), reqsrsapp(is,im), reqsrsapp(is,im)/reqsrs(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', reqis(is,im), reqisapp(is,im), reqisapp(is,im)/reqis(is,im));
    fprintf(' & %8.2e & %8.2e', reqsrs(is,im)/reqis(is,im), reqsrsapp(is,im)/reqisapp(is,im));
    fprintf(' \\\\ \\hline');
    
end
    
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal: RE output: mean
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n mean: Method1 = srs, Method2 = is');
fprintf('\n G_0 = normal(%d,%d): theta(is) = %10.4e', mean0, var0, theta(is));

outrehead;

for im = 1:immax
    
    fprintf('\n %2d & %8.2e', m(is,im), 1-p(is,im));
    fprintf(' & %8.2e & %8.2e', meansum(is,im), meansum(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', remsrs(is,im), remsrs(is,im), remsrs(is,im)/remsrs(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', remis(is,im), remis(is,im), remis(is,im)/remis(is,im));
    fprintf(' & %8.2e & %8.2e', remsrs(is,im)/remis(is,im), remsrs(is,im)/remis(is,im));
    fprintf(' \\\\ \\hline');
    
end
    
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal: RE output: ec
% srs and is
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n ec: Method1 = srs, Method2 = is');
fprintf('\n G_0 = normal(%d,%d): theta(is) = %10.4e', mean0, var0, theta(is));

outrehead;

for im = 1:immax
    
    fprintf('\n %2d & %8.2e', m(is,im), 1-p(is,im));
    fprintf(' & %8.2e & %8.2e', ec(is,im), ecapp(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', reecsrs(is,im), reecsrsapp(is,im), reecsrsapp(is,im)/reecsrs(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', reecis(is,im), reecisapp(is,im), reecisapp(is,im)/reecis(is,im));
    fprintf(' & %8.2e & %8.2e', reecsrs(is,im)/reecis(is,im), reecsrsapp(is,im)/reecisapp(is,im));
    fprintf(' \\\\ \\hline');
    
end
    
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal: RE output: ec
% msis and isdm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n ec: Method1 = msis, Method2 = isdm');
fprintf('\n G_0 = normal(%d,%d): theta(is) = %10.4e', mean0, var0, theta(is));

outrehead;

for im = 1:immax
    
    fprintf('\n %2d & %8.2e', m(is,im), 1-p(is,im));
    fprintf(' & %8.2e & %8.2e', ec(is,im), ecapp(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', reecmsis(is,im), reecmsisapp(is,im), reecmsisapp(is,im)/reecmsis(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', reecisdm(is,im), reecisdmapp(is,im), reecisdmapp(is,im)/reecisdm(is,im));
    fprintf(' & %8.2e & %8.2e', reecsrs(is,im)/reecmsis(is,im), reecsrsapp(is,im)/reecmsisapp(is,im));
    fprintf(' \\\\ \\hline');
    
end
    
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal: RE output: ec
% msis and de
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n ec: Method1 = msis, Method2 = de');
fprintf('\n G_0 = normal(%d,%d): theta(is) = %10.4e', mean0, var0, theta(is));

outrehead;

for im = 1:immax
    
    fprintf('\n %2d & %8.2e', m(is,im), 1-p(is,im));
    fprintf(' & %8.2e & %8.2e', ec(is,im), ecapp(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', reecmsis(is,im), reecmsisapp(is,im), reecmsisapp(is,im)/reecmsis(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', reecde(is,im), reecdeapp(is,im), reecdeapp(is,im)/reecde(is,im));
    fprintf(' & %8.2e & %8.2e', reecsrs(is,im)/reecmsis(is,im), reecsrsapp(is,im)/reecmsisapp(is,im));
    fprintf(' \\\\ \\hline');
    
end
    
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal: RE output table
% srs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n ***********************************************');

fprintf('\n G_0 = normal(%d,%d): theta(is) = %10.4e', mean0, var0, theta(is));

fprintf('\n srs: ec exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reecsrs(is,im));
end
fprintf('\n\n')

fprintf('\n srs: ec app');
fprintf('\n   m        REapp');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reecsrsapp(is,im));
end
fprintf('\n\n')

fprintf('\n srs: q exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reqsrs(is,im));
end
fprintf('\n\n')

fprintf('\n srs: q app');
fprintf('\n   m        REapp');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reqsrsapp(is,im));
end
fprintf('\n\n')

fprintf('\n srs: mean exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), remsrs(is,im));
end
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal: RE output table
% is
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n ***********************************************');
fprintf('\n ***********************************************');

fprintf('\n G_0 = normal(%d,%d): theta(is) = %10.4e', mean0, var0, theta(is));

fprintf('\n is: ec exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reecis(is,im));
end
fprintf('\n\n')

fprintf('\n is: ec app');
fprintf('\n   m        REapp');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reecisapp(is,im));
end
fprintf('\n\n')

fprintf('\n is: q exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reqis(is,im));
end
fprintf('\n\n')

fprintf('\n is: q app');
fprintf('\n   m        REapp');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reqisapp(is,im));
end
fprintf('\n\n')

fprintf('\n is: mean exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), remis(is,im));
end
fprintf('\n\n')

fprintf('\n is: mean exact, var');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), varmis(is,im));
end
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal: RE output table
% msis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n ***********************************************');
fprintf('\n ***********************************************');

fprintf('\n G_0 = normal(%d,%d): theta(is) = %10.4e', mean0, var0, theta(is));

fprintf('\n msis: ec exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reecmsis(is,im));
end
fprintf('\n\n')

fprintf('\n msis: ec app');
fprintf('\n   m        REapp');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reecmsisapp(is,im));
end
fprintf('\n\n')

fprintf('\n msis: q exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reqis(is,im));
end
fprintf('\n\n')

fprintf('\n msis: q app');
fprintf('\n   m        REapp');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reqisapp(is,im));
end
fprintf('\n\n')

fprintf('\n msis: mean exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), remsrs(is,im));
end
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal: RE output table
% isdm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n ***********************************************');
fprintf('\n ***********************************************');

fprintf('\n G_0 = normal(%d,%d): theta(is) = %10.4e', mean0, var0, theta(is));

fprintf('\n isdm: ec exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reecisdm(is,im));
end
fprintf('\n\n')

fprintf('\n isdm: ec app');
fprintf('\n   m        REapp');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reecisdmapp(is,im));
end
fprintf('\n\n')

fprintf('\n isdm: q exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reqisdm(is,im));
end
fprintf('\n\n')

fprintf('\n isdm: q app');
fprintf('\n   m        REapp');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reqisdmapp(is,im));
end
fprintf('\n\n')

fprintf('\n isdm: mean exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), remisdm(is,im));
end
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal: RE output table
% de: double estimator with fixed weights
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n ***********************************************');
fprintf('\n ***********************************************');

fprintf('\n G_0 = normal(%d,%d): theta(is) = %10.4e', mean0, var0, theta(is));

fprintf('\n de: ec exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reecde(is,im));
end
fprintf('\n\n')

fprintf('\n de: ec app');
fprintf('\n   m        REapp');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reecdeapp(is,im));
end
fprintf('\n\n')

fprintf('\n de: q exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reqde(is,im));
end
fprintf('\n\n')

fprintf('\n de: q app');
fprintf('\n   m        REapp');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reqdeapp(is,im));
end
fprintf('\n\n')

fprintf('\n de: mean exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), remde(is,im));
end
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal: RE output table
% deo: double estimator with optimal weights
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n ***********************************************');
fprintf('\n ***********************************************');

fprintf('\n G_0 = normal(%d,%d): theta(is) = %10.4e', mean0, var0, theta(is));

fprintf('\n de optimal: ec exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reecdeo(is,im));
end
fprintf('\n\n')

% fprintf('\n de: ec app');
% fprintf('\n   m        REapp');
% for im = 1:immax
%     fprintf('\n ( %2d , %12.6e )', m(is,im), reecdeapp(is,im));
% end
% fprintf('\n\n')

fprintf('\n deo: q exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), reqdeo(is,im));
end
fprintf('\n\n')

% fprintf('\n de: q app');
% fprintf('\n   m        REapp');
% for im = 1:immax
%     fprintf('\n ( %2d , %12.6e )', m(is,im), reqdeapp(is,im));
% end
% fprintf('\n\n')

fprintf('\n deo: mean exact');
fprintf('\n   m        RE');
for im = 1:immax
    fprintf('\n ( %2d , %12.6e )', m(is,im), remdeo(is,im));
end
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normal
% f(q(is,im)), SD, etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n\n')

fprintf('\n q: G_0 = normal(%d,%d): theta = %12.6e', mean0, var0, theta(is));
fprintf('\n m    1-p');
fprintf('        q          qapp       qapp/q');
fprintf('     f(q)       f(qapp)   fapp(qap)');
fprintf(' fa(qa)/f(qa) fa(qa)/f(q)');
fprintf(' IS:SD     IS:SDapp   SDapp/SD');
fprintf('\n');

fprintf('\n$m$ & $1-p$    ');
fprintf('& $\\q$     & $\\qapp$  & $\\frac{\\qapp}{\\q}$');
fprintf('& $f(\\q)$ & $f(\\qapp)$ & $\\fapp(\\qapp)$ ');
fprintf('& $\\frac{\\fapp(\\qapp)}{f(\\qapp)}$');
fprintf('& $\\frac{\\fapp(\\qapp)}{f(\\q)}$');
fprintf('& $\\sdqis$ & $\\sdqisapp$ & $\\frac{\\sdqisapp}{\\sdqis}$');
fprintf('\n\\\\ \\hline');


% Normal: f(q(is,im)), SD, etc.    

for im = 1:immax
    
    fprintf('\n %2d & %8.2e', m(is,im), 1-p(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', q(is,im), qapp(is,im), qapp(is,im)/q(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', fq(is,im), fqapp(is,im), fappqapp(is,im));
    fprintf(' & %8.2e & %8.2e', fappqapp(is,im)/fqapp(is,im), fappqapp(is,im)/fq(is,im));
    fprintf(' & %8.2e & %8.2e & %8.2e', sdqnum(is,im), sdqnumapp(is,im), sdqnumapp(is,im)/sdqnum(is,im));
    fprintf('\n\\\\ \\hline');
    
%     fprintf('\n m = %d, sdqnum = %e, sdqnumapp = %e', m(is,im), sdqnum(is,im), sdqnumapp(is,im) );
    
end


fprintf('\n\n')
fprintf('\n\n')
 
 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Erlang
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idist = ierlang;
ismax = 6;

% scale is mean of each stage of Erlang
scale = 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Erlang: RE computations

% loop over different number of stages in Erlang
for is = 1:ismax

    % number of stages is power of 2
    s(is) = 2.^(is-1);
    
    % compute mean and variance of iid sum
    mean0 = s(is) .* scale;
    var0 = s(is) .* (scale.^2);
    
    sd0 = sqrt(var0);

% compute the IS twisting parameter using root-finding    
    fun = @(x) -s(is).*scale.*x ./ (1-(scale.*x)) ...
        - s(is).*log(1-(scale.*x)) + beta;
    theta(is) = fzero(fun, 0);

% MGF(theta) for G_0
    mgf0(is) = ( 1 ./ ( 1 - ( scale.*theta(is) ) ) ).^s(is);

% MGF(-theta) for G_0 and derivatives
    mgf0n(is) = (1/( 1 + (scale.*theta(is)) ) ).^s(is);
    mgf0dn(is) = s(is) .* scale .* ((1+(scale.*theta(is))).^(-s(is)-1));
    mgf0d2n(is) = s(is) .* ( s(is) + 1 ) .* (scale.^2) ...
        .* ((1 + (scale.*theta(is)) ).^(-s(is)-2));

% CGF(theta) for G_0 and derivatives
    cgf0(is) = -s(is) .* log( 1 - (scale.*theta(is)) );
    cgf0d(is) = s(is) .* scale ./ (1 - (scale.*theta(is)) );
    cgf0d2(is) = s(is) .* (scale.^2) ./ (( 1 - (scale.*theta(is)) ).^2);
    
% CGF(-theta) for G_0 and derivatives
    cgf0n(is) = -s(is) .* log( 1 + (scale.*theta(is)) );
    cgf0dn(is) = s(is) .* scale ./ ( 1 + (scale.*theta(is)) );
    cgf0d2n(is) = s(is) .* (scale.^2) ./ (( 1 + (scale.*theta(is)) ).^2);
    
    % the new scale parameter for Erlang under IS
    scalenew = scale ./ ( 1 - (scale.*theta(is)) );
%     fprintf('\n scalenew = %e, theta = %e', scalenew, theta(is));
%     fprintf('\n cgf0 = %e, cgf0d = %e, cgf0d2 = %e', ...
%         cgf0(is), cgf0d(is), cgf0d2(is));


% loop through different number of i.i.d. summands
    for im = 1:immax
    
        % number of summands is power of 2
        m(is,im) = 2.^(im-1);
        p(is,im) = 1.0 - exp(-beta.*m(is,im));

        ms = m(is,im).*s(is);

        % compute mean and variance of iid sum
        meansum(is,im) = m(is,im).*mean0;
        varsum(is,im) = m(is,im).*var0;
        sdsum(is,im) = sqrt(varsum(is,im));

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% q: srs exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute true quantile and true density at true quantile
        q(is,im) = gaminv(p(is,im), ms, scale);
        fq(is,im) = gampdf(q(is,im), ms, scale);
        
%         fprintf(' \n is = %d, im = %d, ms = %d, q = %e, fq = %e', ...
%             is, im, ms, q(is,im), fq(is,im) );

% compute exact asymptotic variance of SRS estimator of quantile
        varqsrs(is,im) = p(is,im) .* ( 1 - p(is,im) ) / ( fq(is,im).^2 );
        sdqsrs(is,im) = sqrt( varqsrs(is,im) );
        reqsrs(is,im) = sdqsrs(is,im)/q(is,im);
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% q: srs approx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute approximate quantile, true density at approx quantile,
% and approximate density at approximate quantile
%
        qapp(is,im) = m(is,im).*cgf0d(is);
        fqapp(is,im) = gampdf(qapp(is,im), ms, scale);
        fappqapp(is,im) = (1/sqrt(2.*pi.*m(is,im).*cgf0d2(is))) .* exp(m(is,im).*(cgf0(is) - theta(is).*cgf0d(is)));

        % compute approximate asymptotic variance of SRS quantile estimator
        varqsrsapp(is,im) = p(is,im).*(1-p(is,im))/( fappqapp(is,im).^2 );
        sdqsrsapp(is,im) = sqrt( varqsrsapp(is,im) );
        reqsrsapp(is,im) = sdqsrsapp(is,im)/qapp(is,im);   
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% q: is exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As with the case of normal i.i.d. summands,
% need to be careful in how variances and covariances
% are numerically computed for Erlang summands. 
% See the discussion around line 320 in this file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        z = q(is,im);
        
        a1z = a1erlang(z);
        a2z = a2erlang(z);
        
%         fprintf('\n a1 = %e, a2 = %e', a1z, a2z);
        
        varqnum(is,im) = a1z + a2z;
        sdqnum(is,im) = sqrt(varqnum(is,im));
        
        varqis(is,im) = varqnum(is,im)/( fq(is,im).^2 );
        sdqis(is,im) = sqrt( varqis(is,im) );
        reqis(is,im) = sdqis(is,im)/q(is,im);
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Erlang: RE   
% q: is approximations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        z = qapp(is,im);

        a1z = a1erlang(z);    
        a2z = a2erlang(z);
        
        varqnumapp(is,im) = a1z + a2z;
        sdqnumapp(is,im) = sqrt(varqnumapp(is,im));

        varqisapp(is,im) = varqnumapp(is,im)/( fappqapp(is,im).^2 );
        sdqisapp(is,im) = sqrt( varqisapp(is,im) );

%        varqnum(is,im) = a1z + a2z;
%        sdqnumapp(is,im) = sqrt(varqnum(is,im));
%        sdqisapp(is,im) = sdqnumapp(is,im)/fappqapp(is,im);

        reqisapp(is,im) = sdqisapp(is,im)/qapp(is,im);
    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% mean: srs exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        varmsrs(is,im) = varsum(is,im);
    
        sdmsrs(is,im) = sqrt( varmsrs(is,im) );
        remsrs(is,im) = sdmsrs(is,im)/meansum(is,im);
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% mean: is exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (m(is,im) == 1)    
            varmis(is,im) = mgf0(is).*mgf0d2n(is) - mean0.^2;
        else
            alpha = mgf0(is).*mgf0n(is);
            varmis(is,im) = m(is,im) .* (alpha.^m(is,im)) ...
                .* ( (m(is,im) .* (cgf0dn(is).^2)) + cgf0d2n(is) ) ...
                - ( meansum(is,im).^2 );
        end
        
        sdmis(is,im) = sqrt(varmis(is,im));
        remis(is,im) = sdmis(is,im)/meansum(is,im);
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% ec: srs exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ec(is,im) = q(is,im) - meansum(is,im);
        covecsrs(is,im) = covsrserlang(q(is,im));
        
%         fprintf('\n ec(%d,%d) = %e', s(is), m(im), ec(is,im));
    
        varecsrs(is,im) = varqsrs(is,im) + varmsrs(is,im) ...
            - 2.*covecsrs(is,im)/fq(is,im);
        sdecsrs(is,im) = sqrt( varecsrs(is,im) );
        reecsrs(is,im) = sdecsrs(is,im)/ec(is,im);

        
        corrsrs(is,im) = (covecsrs(is,im)/fq(is,im))/ ...
        sqrt(varqsrs(is,im).*varmsrs(is,im));

%         fprintf('\n Erlang, s = %d, m = %d, corrsrs = %e, corrsrs^2 = %e', ...
%             s(is), m(is,im), corrsrs(is,im),  corrsrs(is,im).^2);
%         fprintf('\n');

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% ec: srs approx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ecapp(is,im) = qapp(is,im) - meansum(is,im);
        covecsrsapp(is,im) = covsrserlang(q(is,im));
    
        varecsrsapp(is,im) = varqsrsapp(is,im) + varmsrs(is,im) ...
            - 2.*covecsrsapp(is,im)/fappqapp(is,im);
        sdecsrsapp(is,im) = sqrt( varecsrsapp(is,im) );
        reecsrsapp(is,im) = sdecsrsapp(is,im)/ecapp(is,im);
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% ec: is exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        covecis(is,im) = coviserlang(q(is,im));
    
        varecis(is,im) = varqis(is,im) + varmis(is,im) ...
            - 2.*covecis(is,im)/fq(is,im);

%         fprintf('\n varecis (%d) = %e, covecis = %e', ...
%             m(is,im), varecis(is,im), covecis(is,im) );
% %        fprintf('\n varecis2(%d) = %10.4e', m(is,im), varecisnormal);
    
        sdecis(is,im) = sqrt(varecis(is,im));
        reecis(is,im) = sdecis(is,im)/ec(is,im);

    
        corris(is,im) = covecis(is,im)/...
        sqrt(varqis(is,im).*varmis(is,im));

%         fprintf('\n Erlang, s = %d, m = %d, corris = %e, corris^2 = %e', ...
%             s(is), m(is,im), corris(is,im),  corris(is,im).^2);
%         fprintf('\n');
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% ec: is approx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        covecisapp(is,im) = coviserlang(qapp(is,im));
    
        varecisapp(is,im) = varqisapp(is,im) + varmis(is,im) ...
            - 2.*covecisapp(is,im)/fappqapp(is,im);
    
        sdecisapp(is,im) = sqrt(varecisapp(is,im));
        reecisapp(is,im) = sdecisapp(is,im)/ecapp(is,im);

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% ec: msis exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        varecmsis(is,im) = varqis(is,im)./msispar ...
            + varmsrs(is,im)./(1-msispar);
        sdecmsis(is,im) = sqrt( varecmsis(is,im) );
        reecmsis(is,im) = sdecmsis(is,im)/ec(is,im);
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% ec: msis approx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        varecmsisapp(is,im) = varqisapp(is,im)/msispar ...
            + varmsrs(is,im)/(1-msispar);
        sdecmsisapp(is,im) = sqrt( varecmsisapp(is,im) );
        reecmsisapp(is,im) = sdecmsisapp(is,im)/ecapp(is,im);
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% ec: isdm exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        z = q(is,im);
        den = fq(is,im);

        varqisdm(is,im) = varqnisdmerlang(z) ./ (den.^2);
        varmisdm(is,im) = varmisdmerlang;
        covisdm(is,im) = covisdmerlang(z);

        varecisdm(is,im) = varqisdm(is,im) + varmisdm(is,im) ...
            - 2.*covisdm(is,im) ./ den;
        sdecisdm(is,im) = sqrt( varecisdm(is,im) );
        reecisdm(is,im) = sdecisdm(is,im)/ec(is,im);

        sdqisdm(is,im) = sqrt( varqisdm(is,im) );
        reqisdm(is,im) = sdqisdm(is,im)/q(is,im);

        sdmisdm(is,im) = sqrt( varmisdm(is,im) );
        remisdm(is,im) = sdmisdm(is,im)/meansum(is,im);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% ec: isdm approx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        z = qapp(is,im);
        den = fappqapp(is,im);

        varqisdmapp(is,im) = varqnisdmerlang(z) ./ (den.^2);
        varmisdmapp(is,im) = varmisdmerlang;
        covisdmapp(is,im) = covisdmerlang(z);

        varecisdmapp(is,im) = varqisdmapp(is,im) + varmisdmapp(is,im) ...
            - 2.*covisdmapp(is,im) ./ den;
        sdecisdmapp(is,im) = sqrt( varecisdmapp(is,im) );
        reecisdmapp(is,im) = sdecisdmapp(is,im)/ecapp(is,im);

        sdqisdmapp(is,im) = sqrt( varqisdmapp(is,im) );
        reqisdmapp(is,im) = sdqisdmapp(is,im)/qapp(is,im);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% ec: de exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        z = q(is,im);
        den = fq(is,im);

        temp1 = (1 ./ msispar ) ...
            .* ( (dewtq.^2) .* varqis(is,im) ...
            + (dewtm.^2) .* varmis(is,im) ...
            - 2 .* dewtq .* dewtm .* covecis(is,im) ./ den );

        temp2 = (1 ./ (1 - msispar) ) ...
            .* ( ((1-dewtq).^2) .* varqsrs(is,im) ...
            + ((1-dewtm).^2) .* varmsrs(is,im) ...
            + 2 .* (1-dewtq) .* (1-dewtm) .* covecsrs(is,im) ./ den );

        varqde(is,im) = (1 ./ msispar ) .* (dewtq.^2) .* varqis(is,im) ...
            + (1 ./ (1 - msispar) ) .* ((1-dewtq).^2) .* varqsrs(is,im);

        varmde(is,im) = (1 ./ msispar ) .* (dewtm.^2) .* varmis(is,im) ...
            + (1 ./ (1 - msispar) ) .* ((1-dewtm).^2) .* varmsrs(is,im);

        covde(is,im) = (1 ./ msispar ).* dewtq .* dewtm .* covecis(is,im) ...
            + (1 ./ (1 - msispar) ) .* (1-dewtq) .* (1-dewtm) .* covecsrs(is,im);

        varecde(is,im) = varqde(is,im) + varmde(is,im) ...
            - 2 .* covde(is,im) / den;
        sdecde(is,im) = sqrt( varecde(is,im) );
        reecde(is,im) = sdecde(is,im)/ec(is,im);

        sdqde(is,im) = sqrt( varqde(is,im) );
        reqde(is,im) = sdqde(is,im)/q(is,im);

        sdmde(is,im) = sqrt( varmde(is,im) );
        remde(is,im) = sdmde(is,im)/meansum(is,im);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE   
% ec: de approx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        z = qapp(is,im);
        den = fappqapp(is,im);

        temp1 = (1 ./ msispar ) ...
            .* ( (dewtq.^2) .* varqisapp(is,im) ...
            + (dewtm.^2) .* varmis(is,im) ...
            - 2 .* dewtq .* dewtm .* covecisapp(is,im) ./ den );

        temp2 = (1 ./ (1 - msispar) ) ...
            .* ( ((1-dewtq).^2) .* varqsrsapp(is,im) ...
            + ((1-dewtm).^2) .* varmsrs(is,im) ...
            - 2 .* (1-dewtq) .* (1-dewtm) .* covecsrsapp(is,im) ./ den );

        varqdeapp(is,im) = (1 ./ msispar ) .* (dewtq.^2) .* varqisapp(is,im) ...
            + (1 ./ (1 - msispar) ) .* ((1-dewtq).^2) .* varqsrsapp(is,im);

        varmdeapp(is,im) = (1 ./ msispar ) .* (dewtm.^2) .* varmis(is,im) ...
            + (1 ./ (1 - msispar) ) .* ((1-dewtm).^2) .* varmsrs(is,im);

        covdeapp(is,im) = (1 ./ msispar ).* dewtq .* dewtm .* covecisapp(is,im) ...
            + (1 ./ (1 - msispar) ) .* (1-dewtq) .* (1-dewtm) .* covecsrsapp(is,im);

        varecdeapp(is,im) = varqdeapp(is,im) + varmdeapp(is,im) ...
            - 2 .* covdeapp(is,im) ./ den;
        sdecdeapp(is,im) = sqrt( varecdeapp(is,im) );
        reecdeapp(is,im) = sdecdeapp(is,im)/ecapp(is,im);

        sdqdeapp(is,im) = sqrt( varqdeapp(is,im) );
        reqdeapp(is,im) = sdqdeapp(is,im)/qapp(is,im);
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE
% ec: de optimal exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        z = q(is,im);
        den = fq(is,im);
        
        vqis = varqis(is,im) ./ (msispar);
        vqsrs = varqsrs(is,im) ./ ((1 - msispar));
        vmis = varmis(is,im) ./ msispar;
        vmsrs = varmsrs(is,im) ./ (1 - msispar);
        cis = covecis(is,im) ./ (msispar .* den);
        csrs = covecsrs(is,im) ./ ((1 - msispar) .* den);
        
        dea0 = vqsrs.*vmis - cis.*cis - 2.*cis.*csrs - csrs.*csrs ...
            + vqis.*vmis + vqis.*vmsrs + vqsrs.*vmsrs;
        
        dea1 = vqsrs.*vmis + vqsrs.*vmsrs - vmis.*csrs + vmsrs.*cis ...
            - cis.*csrs - csrs.*csrs;
        
        dea2 = vqis.*vmsrs + vqsrs.*vmsrs - vqis.*csrs + vqsrs.*cis ...
            - cis.*csrs - csrs.*csrs;
        
        dewtqo = dea1 ./ dea0;
        dewtmo = dea2 ./ dea0;
        
        %     temp1 = (1 ./ msispar ) ...
        %         .* ( (dewtqo.^2) .* varqis(is,im) ...
        %         + (dewtmo.^2) .* varmis(is,im) ...
        %         - 2 .* dewtqo .* dewtmo .* covecis(is,im) ./ den );
        %
        %     temp2 = (1 ./ (1 - msispar) ) ...
        %         .* ( ((1-dewtqo).^2) .* varqsrs(is,im) ...
        %         + ((1-dewtmo).^2) .* varmsrs(is,im) ...
        %         - 2 .* (1-dewtqo) .* (1-dewtmo) .* covecsrs(is,im) ./ den );
        
        varqdeo(is,im) = (1 ./ msispar ) .* (dewtqo.^2) .* varqis(is,im) ...
            + (1 ./ (1 - msispar) ) .* ((1-dewtqo).^2) .* varqsrs(is,im);
        
        varmdeo(is,im) = (1 ./ msispar ) .* (dewtmo.^2) .* varmis(is,im) ...
            + (1 ./ (1 - msispar) ) .* ((1-dewtmo).^2) .* varmsrs(is,im);
        
        covdeo(is,im) = (1 ./ msispar ).* dewtqo .* dewtmo .* covecis(is,im) ...
            + (1 ./ (1 - msispar) ) .* (1-dewtqo) .* (1-dewtmo) .* covecsrs(is,im);
        
        varecdeo(is,im) = varqdeo(is,im) + varmdeo(is,im) ...
            - 2 .* covdeo(is,im) / den;
        sdecdeo(is,im) = sqrt( varecdeo(is,im) );
        reecdeo(is,im) = sdecdeo(is,im)/abs(ec(is,im));
        
        sdqdeo(is,im) = sqrt( varqdeo(is,im) );
        reqdeo(is,im) = sdqdeo(is,im)/abs(q(is,im));
        
        sdmdeo(is,im) = sqrt( varmdeo(is,im) );
        remdeo(is,im) = sdmdeo(is,im)/abs(meansum(is,im));
        
        
        
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erlang: RE output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n\n')

fprintf('\n G_0 = Erlang');

outreheadq;

for is = 1:ismax

    fprintf('\n\\hline');
    fprintf('\n \\multicolumn{12}{|c|}');
    fprintf('{Erlang:  \\quad $s = %2d$, \\quad $\\tpq =$ %12.6e}', s(is), theta(is));
    fprintf('\n\\\\ \\hline');

    for im = 1:immax
    
        fprintf('\n %2d & %8.2e', m(is,im), 1-p(is,im));
        fprintf(' & %8.2e & %8.2e', q(is,im), qapp(is,im));
        fprintf(' & %8.2e & %8.2e & %8.2e', reqsrs(is,im), reqsrsapp(is,im), reqsrsapp(is,im)/reqsrs(is,im));
        fprintf(' & %8.2e & %8.2e & %8.2e', reqis(is,im), reqisapp(is,im), reqisapp(is,im)/reqis(is,im));
        fprintf(' & %8.2e & %8.2e', reqsrs(is,im)/reqis(is,im), reqsrsapp(is,im)/reqisapp(is,im));
        fprintf(' \\\\ \\hline');
        
    end
end

fprintf('\n\n')







fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n\n')

fprintf('\n G_0 = Erlang');

for is = 1:ismax

    fprintf('\n\\hline');
    fprintf('\n \\multicolumn{12}{|c|}');
    fprintf('{Erlang:  \\quad $s = %2d$, \\quad $\\tpq =$ %12.6e}', s(is), theta(is));
    fprintf('\n\\\\ \\hline');

    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Erlang: RE output table
% srs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n ***********************************************');
    fprintf('\n ***********************************************');
    fprintf('\n ***********************************************');

    fprintf('\n G_0 = Erlang(%d,%d): theta(is) = %10.4e', s(is), scale, theta(is));

    fprintf('\n srs: ec exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reecsrs(is,im));
    end
    fprintf('\n\n')

    fprintf('\n srs: ec app');
    fprintf('\n   m        REapp');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reecsrsapp(is,im));
    end
    fprintf('\n\n')

    fprintf('\n srs: q exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reqsrs(is,im));
    end
    fprintf('\n\n')

    fprintf('\n srs: q app');
    fprintf('\n   m        REapp');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reqsrsapp(is,im));
    end
    fprintf('\n\n')

    fprintf('\n srs: mean exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), remsrs(is,im));
    end
    fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Erlang: RE output table
% is
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n ***********************************************');
    fprintf('\n ***********************************************');

    fprintf('\n G_0 = Erlang(%d,%d): theta(is) = %10.4e', s(is), scale, theta(is));

    fprintf('\n is: ec exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reecis(is,im));
    end
    fprintf('\n\n')

    fprintf('\n is: ec app');
    fprintf('\n   m        REapp');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reecisapp(is,im));
    end
    fprintf('\n\n')

    fprintf('\n is: q exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reqis(is,im));
    end
    fprintf('\n\n')

    fprintf('\n is: q app');
    fprintf('\n   m        REapp');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reqisapp(is,im));
    end
    fprintf('\n\n')

    fprintf('\n is: mean exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), remis(is,im));
    end
    fprintf('\n\n')

    fprintf('\n is: mean exact, var');
    fprintf('\n   m        var');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), varmis(is,im));
    end
    fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Erlang: RE output table
% msis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n ***********************************************');
    fprintf('\n ***********************************************');

    fprintf('\n G_0 = Erlang(%d,%d): theta(is) = %14.8e', s(is), scale, theta(is));

    fprintf('\n msis: ec exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reecmsis(is,im));
    end
    fprintf('\n\n')

    fprintf('\n msis: ec app');
    fprintf('\n   m        REapp');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reecmsisapp(is,im));
    end
    fprintf('\n\n')

    fprintf('\n msis: q exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reqis(is,im));
    end
    fprintf('\n\n')

    fprintf('\n msis: q app');
    fprintf('\n   m        REapp');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reqisapp(is,im));
    end
    fprintf('\n\n')

    fprintf('\n msis: mean exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), remsrs(is,im));
    end
    fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Erlang: RE output table
% isdm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n ***********************************************');
    fprintf('\n ***********************************************');

    fprintf('\n G_0 = Erlang(%d,%d): theta(is) = %10.4e', s(is), scale, theta(is));

    fprintf('\n isdm: ec exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reecisdm(is,im));
    end
    fprintf('\n\n')

    fprintf('\n isdm: ec app');
    fprintf('\n   m        REapp');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reecisdmapp(is,im));
    end
    fprintf('\n\n')

    fprintf('\n isdm: q exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reqisdm(is,im));
    end
    fprintf('\n\n')

    fprintf('\n isdm: q app');
    fprintf('\n   m        REapp');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reqisdmapp(is,im));
    end
    fprintf('\n\n')

    fprintf('\n isdm: mean exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), remisdm(is,im));
    end
    fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Erlang: RE output table
% de
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n ***********************************************');
    fprintf('\n ***********************************************');

    fprintf('\n G_0 = Erlang(%d,%d): theta(is) = %10.4e', s(is), scale, theta(is));

    fprintf('\n de: ec exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reecde(is,im));
    end
    fprintf('\n\n')

    fprintf('\n de: ec app');
    fprintf('\n   m        REapp');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reecdeapp(is,im));
    end
    fprintf('\n\n')

    fprintf('\n de: q exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reqde(is,im));
    end
    fprintf('\n\n')

    fprintf('\n de: q app');
    fprintf('\n   m        REapp');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reqdeapp(is,im));
    end
    fprintf('\n\n')

    fprintf('\n de: mean exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), remde(is,im));
    end
    fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Erlang: RE output table
% de optimal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n ***********************************************');
    fprintf('\n ***********************************************');

    fprintf('\n G_0 = Erlang(%d,%d): theta(is) = %10.4e', s(is), scale, theta(is));

    fprintf('\n deo: ec exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reecdeo(is,im));
    end
    fprintf('\n\n')

%     fprintf('\n de: ec app');
%     fprintf('\n   m        REapp');
%     for im = 1:immax
%         fprintf('\n ( %2d , %12.6e )', m(is,im), reecdeapp(is,im));
%     end
%     fprintf('\n\n')

    fprintf('\n deo: q exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), reqdeo(is,im));
    end
    fprintf('\n\n')

%     fprintf('\n de: q app');
%     fprintf('\n   m        REapp');
%     for im = 1:immax
%         fprintf('\n ( %2d , %12.6e )', m(is,im), reqdeapp(is,im));
%     end
%     fprintf('\n\n')

    fprintf('\n deo: mean exact');
    fprintf('\n   m        RE');
    for im = 1:immax
        fprintf('\n ( %2d , %12.6e )', m(is,im), remdeo(is,im));
    end
    fprintf('\n\n')


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Erlang:
% f(q(is,im)), SD, etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n ***********************************************');
fprintf('\n ***********************************************');
fprintf('\n\n')

fprintf('\n G_0 = Erlang');
fprintf('\n m    1-p');
fprintf('        q          qapp       qapp/q');
fprintf('     f(q)       f(qapp)   fapp(qap)');
fprintf(' fa(qa)/f(qa) fa(qa)/f(q)');
fprintf(' IS:SD     IS:SDapp   SDapp/SD');
fprintf('\n');

fprintf('\n$m$ & $1-p$    ');
fprintf('& $\\q$     & $\\qapp$  & $\\frac{\\qapp}{\\q}$');
fprintf('& $f(\\q)$ & $f(\\qapp)$ & $\\fapp(\\qapp)$ ');
fprintf('& $\\frac{\\fapp(\\qapp)}{f(\\qapp)}$');
fprintf('& $\\frac{\\fapp(\\qapp)}{f(\\q)}$');
fprintf('& $\\sdqis$ & $\\sdqisapp$ & $\\frac{\\sdqisapp}{\\sdqis}$');
fprintf('\n\\\\ \\hline');




% scale = 1;

for is = 1:6

    fprintf('\n\\hline');
    fprintf('\n \\multicolumn{13}{|c|}');
    fprintf('{Erlang:  \\quad $s = %2d$, \\quad $\\tpq =$ %12.6e}', s(is), theta(is));
    fprintf('\n\\\\ \\hline');

    
    for im = 1:immax
    
        fprintf('\n %2d & %8.2e', m(is,im), 1-p(is,im));
        fprintf(' & %8.2e & %8.2e & %8.2e', q(is,im), qapp(is,im), qapp(is,im)/q(is,im));
        fprintf(' & %8.2e & %8.2e & %8.2e', fq(is,im), fqapp(is,im), fappqapp(is,im));
        fprintf(' & %8.2e & %8.2e', fappqapp(is,im)/fqapp(is,im), fappqapp(is,im)/fq(is,im));
        fprintf(' & %8.2e & %8.2e & %8.2e', sdqnum(is,im), sdqnumapp(is,im), sdqnumapp(is,im)/sdqnum(is,im));
        fprintf('\n\\\\ \\hline');

    
    end
    
end


fprintf('\n\n')
 



end


