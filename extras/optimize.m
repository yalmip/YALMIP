function varargout = optimize(varargin)
%OPTIMIZE Computes solution to optimization problem
%
%   DIAGNOSTIC = OPTIMIZE(Constraint,Objective,options) is the common command to
%   solve optimization problems of the following kind
%
%    min        Objective
%    subject to
%            Constraint
%
%   NOTES
%    To obtain solution for a variable, use VALUE.
%    To obtain dual variable for a constraint, use DUAL.
%    See YALMIPERROR for error codes returned in output.
%
%   OUTPUT
%     diagnostic : Diagnostic information
%
%   INPUT
%     Constraint : Object describing constraints. Can be [].
%     Objective  : Object describing the objective. Can be [].
%     options    : Options structure. See SDPSETTINGS. Can be [].
%
%   EXAMPLE
%    A = randn(15,5);b = rand(15,1)*5;c = randn(5,1);
%    x = sdpvar(5,1);
%    diagnostics = optimize(A*x<=b,c'*x);value(x)
%
%   See also SDPVAR, @SDPVAR/VALUE, SDPSETTINGS, YALMIPERROR

[varargout{1:nargout}] = solvesdp(varargin{:});
