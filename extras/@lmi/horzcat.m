function F = vertcat(varargin)

F = set(varargin{1});
for i=2:1:nargin
    F=F+set(varargin{i});
end