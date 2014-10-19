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
%
% If you intend to use the hull operator in a parameterized fashion in an 
% optimizer setting, you can declare the parameters in the last argument.
% These variables will then be considered as if the were constants
%
% H = hull(F1,F2,...,Parameters)



% Pull out the parameter if one is given
ParameterVariables = [];
if isa(varargin{end},'sdpvar')
    Parameter = varargin{end};
    ParameterVariables = getvariables(Parameter);
    varargin = {varargin{1:end-1}};     
end
N = length(varargin);

if N==1
    Fhull = varargin{1};
    t = [];
    y = [];
end

% Pre-process to convert convex quadratic constraints to socp constraints.
% This makes the perspective code easier.
% We also wexpand all graph-operators etc to find out all involved
% variables and constraints
for i = 1:N
    varargin{i} = convertquadratics(varargin{i});
    varargin{i} = expandmodel(varargin{i},[]);
end

variables = [];
for i = 1:N
    % quick fix
    if isa(varargin{i},'constraint')
        varargin{i} = lmi(varargin{i});
    end
    if ~(isa(varargin{i},'lmi') | isa(varargin{i},'socc'))
        error('Hull can only be applied to linear constraints');
    elseif ~(islinear(varargin{i}))
        error('Hull can only be applied to linear and convex quadratic constraints');
    end
    variables = unique([variables depends(varargin{i})]);
    variables = setdiff(variables,ParameterVariables);
end

if N == 1
    Fhull = varargin{1};
    return
end

% Define variables. To allow users to easily plot convex hulls, we mark the
% internally defined variables as auxilliary. YALMIP will then not plot
% w.r.t these variables
nInitial = yalmip('nvars');
y = sdpvar(repmat(length(variables),1,N),repmat(1,1,N));
t = sdpvar(N,1);
nNow = yalmip('nvars');
yalmip('addauxvariables',nInitial+1:nNow);

Fhull = lmi([]);
for i = 1:N
    Fi = flatten(varargin{i});
    tvariable = getvariables(t(i));
    for j = 1:length(Fi.clauses)       
        Xi = Fi.clauses{j}.data;
        if ~isempty(ParameterVariables)
            Xitrue = Xi;
            Xi = replace(Xi,Parameter,0);
        end
        local_variables = getvariables(Xi);
        local_index = find(ismember(variables,local_variables));
        new_variables = getvariables(y{i}(local_index));
        Xipersp = brutepersp(Xi,tvariable,new_variables); 
        if ~isempty(ParameterVariables)
            Fi.clauses{j}.data = Xipersp + t(i)*(Xitrue-Xi);
        else
            Fi.clauses{j}.data = Xipersp;
        end
        Fi.clauses{j}.handle = ['F(y_' num2str(i) ')'];
    end
    Fhull = Fhull + Fi;
end
Fhull = Fhull + [(sum([y{:}],2) == recover(variables)):'sum y_i == x'];
Fhull = Fhull + [(sum(t)==1):'Multipliers sum to 1'] + [(t>=0):'Positive multiplier'];
Fhull = expanded(Fhull,1);
yalmip('setdependence',[reshape([y{:}],[],1);t(:)],recover(variables));
%yalmip('setdependence',recover([10 11]),recover(variables));
yalmip('addauxvariables',getvariables([reshape([y{:}],[],1);t(:)]));
y = [y{:}];