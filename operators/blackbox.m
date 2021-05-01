function y = blackbox(varargin)
%BLACKBOX Define user-specified function
%
% y = blackbox(f,x,df)
%
% f     : Function handle 
% x     : SDPVAR
% df    : Optional functional handle
%
% See also SDPVAR

% This is just rebranding...

if nargin == 2
    y = sdpfun(varargin{2},varargin{1});
elseif nargin == 3
    y = sdpfun(varargin{2},varargin{1},varargin{3});
else
    error('The number of arguments to BLACKBOX should be 2 or 3 (f,x,df)')
end
