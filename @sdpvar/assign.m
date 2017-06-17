function assign(X,value,ls)
%ASSIGN Assigns a numerical value to an sdpvar
%
% ASSIGN(X,value)   Tries to set the free variables in X so that
%                   double(X)=value. Notice that other variables
%                   sharing the same free variables will be affected.
%                   If the assignment is infeasible, an error message
%                   will be issued.
%
% ASSIGN(X,value,1) Least square assignment.
%
% After assigning values to decision variables, these can be communicated
% to the solver as an initial guess by turning on the option 'usex0' in the
% options structure sent to optimize.

if nargin<3
    ls = 0;
end

if ~isa(X,'sdpvar')
    error('First argument should be an SDPVAR object.');
end

if isa(value,'logical')
    value = double(value);
end

if ~isnumeric(value)
    error('Second argument should be NUMERIC.');
end

if prod(size(value)) == 1
    value = repmat(value,size(X));
end

if isempty(value)
    return
end

if ~isequal(size(X),size(value))
    error('Both arguments must have same size')
end
if ~isa(X,'sdpvar')
    error('First arguments must be an sdpvar object')
end

if is(X,'complex')
    assign(real(X),real(value));
    assign(imag(X),imag(value));
    return
end

x_lmi_variables = X.lmi_variables;
b = value(:)-X.basis(:,1);
A = X.basis(:,2:end);
feas_var = A\b;
% Improve
e = A*feas_var-b;
de = A\e;
feas_var = feas_var-de;

if ~ls
    if norm(A*feas_var-b)>sqrt(eps)
        error('Inconsistent assignment')
    end
end

sol = yalmip('getsolution');
keep_these = find(~ismember(sol.variables,x_lmi_variables));
sol.optvar = [sol.optvar(keep_these);feas_var(:)];
sol.variables = [sol.variables(keep_these);x_lmi_variables(:)];
yalmip('setallsolution',sol);



