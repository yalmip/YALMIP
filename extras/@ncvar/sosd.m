function [h,parvar,Q,v] = sosd(X)
%SOSD Returns sum-of-squares decomposition
%
% [h,parvar,Q,v] = sosd(X)
%
% Example
%
%  x = sdpvar(1,1);F = (sos(1+x^2+x^4));solvesos(F);
%  sdisplay(sosd(F))

[h,parvar,Q,v] = yalmip('getsos',X.extra.sosid);    

