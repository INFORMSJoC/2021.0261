%
% function to compute likelihood ratio for ISDM.
% 
%
function lr = lrisdm(y)

global m im theta beta
global s is ms scale
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global msispar dewtq dewtm
global p


temp = msispar .* exp( theta(is).*y - m(is,im).*cgf0(is) ) + (1-msispar);
lr = 1 ./ temp;


end


