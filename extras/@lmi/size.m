function varargout=size(varargin)
%size              Returns the number of constraints
%
%    n = SIZE(F)     Returns the number of constraints
%    [n,m] = SIZE(F) Returns the number of constraints, m=1

F = flatten(varargin{1});
switch (nargout)
    case {0,1}
        varargout{1} = [length(F.clauses) 1];
    case 2
        varargout{1} = length(F.clauses);
        varargout{2} = 1;
    otherwise
        error('>2 outputs in size?');
end
