function [h,parvar,Q,v] = sosd(X,ops)
%SOSD Returns sum-of-squares decomposition
%
% h = sosd(p)
%
% Example
%
%  x = sdpvar(1,1);sdisplay(sosd(x^4 + 1))

if isfield(X.extra,'sosid')
    % This is used internally when calling sosd on SOS constraint
    [h,parvar,Q,v] = yalmip('getsos',X.extra.sosid);
else
    F = sos(X);
    if nargin == 1
        solvesos(F,[],sdpsettings('verbose',0));
    else
        solvesos(F,[],ops);
    end
    [h,parvar,Q,v] = sosd(F);
end
