function N=numel(varargin)
%NUMEL (overloaded)

% Author Johan Löfberg 
% $Id: numel.m,v 1.4 2006-07-26 20:17:58 joloef Exp $   

X = varargin{1};
N = prod(size(X));
