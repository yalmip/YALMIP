function y=any(varargin)
%ANY (overloaded)

% Author Johan Löfberg 
% $Id: any.m,v 1.1 2005-10-12 16:05:54 joloef Exp $   

if nargin == 1
    y = any(varargin{1}.basis);
else
    X = varargin{1};
    p = varargin{2};
    y = any(X.basis,p);
end
