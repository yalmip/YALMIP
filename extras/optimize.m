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
%    To obtain solution for a variable, use DOUBLE.
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
%    diagnostics = optimize(A*x<=b,c'*x);double(x)
%
%   See also SDPVAR, @SDPVAR/DOUBLE, SDPSETTINGS, YALMIPERROR

% Author Johan Löfberg
% $Id: optimize.m,v 1.75 2010-04-06 06:32:49 joloef Exp $

varargout{:} = solvesdp(varargin{:});
