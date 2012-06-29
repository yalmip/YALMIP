function varargout = find(X)
%FIND (overloaded)

% Author Johan Löfberg
% $Id: find.m,v 1.1 2005-10-12 16:05:54 joloef Exp $

[varargout{1:nargout}] = find(X.basis);