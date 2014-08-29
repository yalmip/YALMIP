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

if nargin<3
    ls = 0;
end

if ~isa(X,'ncvar')
    error('First argument should be an NCVAR object.');
end

if ~isa(value,'double')
    error('Second argument should be a DOUBLE.');
end

if prod(size(value)) == 1
    value = repmat(value,size(X));
end

if ~isequal(size(X),size(value))
  error('Both arguments must have same size') 
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



