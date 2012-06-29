function X = ndsdpvar(varargin)
% NDSDPVAR Constructor for multi-dimensional SDPVAR object

% Author Johan Löfberg
% $Id: ndsdpvar.m,v 1.10 2007-01-04 08:52:49 joloef Exp $

% Sometimes it is convenient to cast a 2D SDPVAR variables as an nD
% variable

if nargin == 1 & isa(varargin{1},'sdpvar')
    X = varargin{1};
    X = struct(X);
    X.conicinfo = [0 0];
    X = class(X,'ndsdpvar');
    return
end

type  = 'symmetric';
field = 'real';

d = 0;
i = 1;
while i<=nargin
    if isa(varargin{i},'double')
        d = d + 1;
    end
    i = i + 1;
end

n = [varargin{1:d}];

if nargin > d
    type = varargin{d+1};
    if n(1)~=n(2) & ~isempty(strmatch(type,'symmetric'))
        error('non-square matrix cannot be symmetric');
    end
else
    if n(1)==n(2)
        type = 'symmetric';
    else
        type = 'full';
    end
end
if nargin > d+1
    field = varargin{d+2};
end

X = [];
for i = 1:prod(n(3:end))
    x = sdpvar(n(1),n(2),type,field);
    X = [X;x(:)];
end
X = struct(X);
X.dim = n;
X = class(X,'ndsdpvar');
X = clean(X);