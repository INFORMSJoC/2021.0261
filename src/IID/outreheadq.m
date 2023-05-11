%
% function to numerically compute first term in
% exact and approximate variance of importance-sampling 
% estimator of p-quantile of sum of m i.i.d. random variables 
% with marginal CDF G_0 as N(0,1).
% 
%
function outreheadq

global m im theta beta
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global p
fprintf('\n                                ');
fprintf('       SRS        SRS        SRS       ');
fprintf(' IS         IS         IS');
fprintf('         RE         REapp');

fprintf('\n m    1-p        measure    app ');
fprintf('       RE         REapp      REapp/RE  ');
fprintf(' RE         REapp      REapp/RE');
fprintf('   SRS/IS     SRS/IS');
fprintf('\n');

fprintf('\n$m$ & $1-p$    & $\\q$     & $\\qapp$');
fprintf(' & $\\REsrs$ & $\\REsrsapp$ & $\\frac{\\REsrsapp}{\\REsrs}$');
fprintf('& $\\REis$ & $\\REisapp$ & $\\frac{\\REisapp}{\\REis}$');
fprintf('& $\\frac{\\REsrs}{\\REis}$ & $\\frac{\\REsrsapp}{\\REisapp}$');
fprintf('\n\\\\ \\hline');

end