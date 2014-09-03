function sys = optproblem(Constraints,Objective,Options)
%OPTPROBLEM  Container for optimization problem
%
%   P = OPTPROBLEM(Constraints,Objective,Options) collects constraints,
%   objective and options of an optimization problem
%
%   Example
%
%    The following code creates an optimization problem, and then minimizes
%    the objective function 
%
%    x = sdpvar(1);P = optproblem(x >= 0, x^2);minimize(P)

if nargin < 2
    Objective = [];
end
if nargin < 3
    Options = sdpsettings;
end

if isa(Constraints,'constraint')
    Constraints = lmi(Constraints);
end
sys.Constraints = Constraints;
sys.Objective = Objective;
sys.Options = Options;
sys = class(sys,'optproblem');

