%
% function to compute likelihood ratio for IS
% (not ISDM).
% 
%
function out = lr(y)

global m im theta beta
global s is ms scale
global meansum sdsum
global mgf0 cgf0 cgf0d cgf0d2 q qapp
global msispar dewtq dewtm
global p


out = exp( m(is,im).*cgf0(is) - theta(is).*y );

end


