function [Fhull,t,y] = hull(varargin)
% HULL  Construct a model of the convex hull
%
% H = hull(F1,F2,...)
%
% OUTPUT
%   H   : Constraint object describing the convex hull of the input constraints
%
% INPUT
%   Fi  : Constraint objects with constraints
%
% Note that the convex representation of the convex hull requires a lifting
% (introduction of auxially variables). Hence, if you have many set of
% constraints, your problem rapidly grows large.

if nargin==1
    Fhull = varargin{1};
    t = [];
end
% Pre-process to convert convex quadratic constraints to socp constraints.
% This makes the perspective code easier.
% We also wexpand all graph-operators etc to find out all involved
% variables and constraints

for i = 1:nargin
    varargin{i} = convertquadratics(varargin{i});
    varargin{i} = expandmodel(varargin{i},[]);
end

variables = [];
for i = 1:nargin
    % quick fix
    if isa(varargin{i},'constraint')
        varargin{i} = set(varargin{i});
    end
    if ~(isa(varargin{i},'lmi') | isa(varargin{i},'socc'))
        error('Hull can only be applied to linear constraints');
    elseif ~(islinear(varargin{i}))
        error('Hull can only be applied to linear and convex quadratic constraints');
    end
    variables = unique([variables depends(varargin{i})]);
end

if nargin == 1
    Fhull = varargin{1};
    return
end

% Define variables. To allow users to easily plot convex hulls, we mark the
% internally defined variables as auxilliary. YALMIP will then not plot
% w.r.t these variables
nInitial = yalmip('nvars');
y = sdpvar(repmat(length(variables),1,nargin),repmat(1,1,nargin));
t = sdpvar(nargin,1);
nNow = yalmip('nvars');
yalmip('addauxvariables',nInitial+1:nNow);

Fhull = set([]);
for i = 1:nargin
    Fi = varargin{i};
    tvariable = getvariables(t(i));
    for j = 1:length(Fi.clauses)
        local_variables = getvariables(Fi);
        Xi = Fi.clauses{j}.data;
        local_variables = getvariables(Xi);
        local_index = find(ismember(variables,local_variables));
        new_variables = getvariables(y{i}(local_index));
        Fi.clauses{j}.data = brutepersp(Fi.clauses{j}.data,tvariable,new_variables);       
    end
    Fhull = Fhull + Fi;
end
Fhull = Fhull + set(sum([y{:}],2) == recover(variables));
Fhull = Fhull + set(sum(t)==1) + set(t>=0);
Fhull = expanded(Fhull,1);
yalmip('setdependence',[reshape([y{:}],[],1);t(:)],recover(variables));
%yalmip('setdependence',recover([10 11]),recover(variables));
yalmip('addauxvariables',getvariables([reshape([y{:}],[],1);t(:)]));
y = [y{:}];