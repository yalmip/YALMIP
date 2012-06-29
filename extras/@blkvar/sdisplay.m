function F = sdisplay(X)
% DISPLAY Overloaded

% Author Johan Löfberg
% $Id: sdisplay.m,v 1.1 2008-05-08 13:41:06 joloef Exp $

try
    error('Not implemented. Convert to SDPVAR first')
%    X = sdpvar(X);
%    evalin('base',sdisplay(X));
catch
    disp('Incomplete block variable.');
end
