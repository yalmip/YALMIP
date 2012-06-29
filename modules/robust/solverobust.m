function sol = solverobust(varargin)
%SOLVEROBUST Solve robust conic problem
%
%    sol = solverobust(F,h,options,uncertain) is used for finding a robust
%    solution to an uncertain conic problem.

%    There is typically no need to use SOLVEROBUST. Use solvesdp directly
%    and simply declare some variables uncertain using UNCERTAIN
%
%    Currently, YALMIP supports robust solutions to conic problems where
%    uncertainty only enter in linear inequality constraints. The
%    uncertainty set is allowed to be an arbitrary conic set. For uncertain
%    SOCP and SDP constraints, the uncertainty set has to be polytopic. 
%
%    The robust problem supported can be formulated as
%
%              min max_w f(x)+c(w)'*x
%
%       subject to  H(x) >(=) 0
%                   A(w)x <= b(w)  for all w:G(w) >(=) 0
%
%    The data c, A and b are linearly parameterized by uncertainty w,
%    the constraints H and G are general conic sets.
%
%   INPUT
%    F         : Object with constraints and uncertainty description
%    h         : scalar SDPVAR object (can be [])
%    options   : options structure obtained from SDPSETTINGS (can be [])
%    uncertain : SDPVAR object defining uncertain variables
%
%   OUTPUT
%    sol       : Solution diagnostic.
%
%   EXAMPLE
%    sdpvar x w
%    F = [x + w <= 1, -0.5 <= w <= 0.5];
%    solverobust(F,-x,[],w) % Optimal value x=0.5
%
%    Recommended version
%    sdpvar x w
%    F = [x + w <= 1, -0.5 <= w <= 0.5, uncertain(w)];
%    solvesdp(F,-x)
%
%   NOTES
%
%    The constraints and objective have to satisfy a number of conditions
%    for  the robustification to be tractable. Please refer to the YALMIP
%    Wiki for the current assumptions (this is constantly developing)
%
% See also SOLVESDP, SOLVESOS, SDPSETTINGS, SDPVAR, SET

if nargin < 3
    varargin{3} = sdpsettings;
    ops = varargin{3};
else
    ops = varargin{3};
    if isempty(ops)
        ops = sdpsettings;
        varargin{3} = ops;
    end
end

% convert to robust model
[F,h,failure] = robustify(varargin{:});
if failure
    error('Failed to create robustified model. Check the Wiki!')
else
    % and solve this model instead
    sol = solvesdp(F,h,ops);
end
