function P = horzcat(varargin)
%MINIMIZSE  Adds constraint to optimization problem
%
%   P = [P,Constraint]


P = varargin{1};
for i = 2:nargin
    if isa(varargin{i},'optproblem')
        P.Constraints = [P.Constraints, varargin{i}.Constraints];
    else
        P.Constraints = [P.Constraints, varargin{i}];
    end
end
