function setsdpvar(X,value,ls)
%SETSDPVAR Obsolete, use ASSIGN instead

% Author Johan Löfberg 
% $Id: setsdpvar.m,v 1.1 2006-08-10 18:00:22 joloef Exp $  

if nargin<3
    ls = 0;
end

if ~isa(X,'sdpvar')
    error('First argument should be an SDPVAR object.');
end

if ~isa(value,'double')
    error('Second argument should be a DOUBLE.');
end

if ~isequal(size(X),size(value))
  error('Both arguments must have same size') 
end
if ~isa(X,'sdpvar')
  error('First arguments must be an sdpvar object') 
end

x_lmi_variables = X.lmi_variables;
b = value(:)-X.basis(:,1);
A = X.basis(:,2:end);
feas_var = A\b;

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



