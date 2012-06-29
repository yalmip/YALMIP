function [h,parvar,Q,v] = sosd(X)
%SOSD Returns sum-of-squares decomposition
%
% [h,parvar,Q,v] = sosd(X)
%
% Example
%
%  x = sdpvar(1,1);F = set(sos(1+x^2+x^4));solvesos(F);
%  sdisplay(sosd(F))

%if ~is(X,'sos')
%    error('SOS decompositions are only available for SOS objects');
%else
[h,parvar,Q,v] = yalmip('getsos',X.extra.sosid);    
%end
