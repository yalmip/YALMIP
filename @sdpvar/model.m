function [properties,F,arguments,fcn]=model(X,method,options,extstruct,w)
%MODEL  Internal function to extracts nonlinear operator models
%
% [properties,F] = model(x)
%
% MODEL returns the constraints needed to model a variable related to an
% extended operator such as min, max, abs, norm, geomean, ...
%
% Examples :
%
% sdpvar x y;
% t = min(x,y);
% [properties,F] = model(t)
% Gives F = [t<=x, t<=y]

extvar = getvariables(X);
arguments   = cell(1,length(extvar));
properties  = cell(1,length(extvar));

if nargin<2
    method = 'graph';
end

if nargin < 3
    options = [];
end

if nargin < 5
    % Not used
    w = [];
end

if nargin<4
    extstruct = yalmip('extstruct',extvar);
elseif isempty(extstruct)
    extstruct = yalmip('extstruct',extvar);
end

if isempty(extstruct)
    error('This is not a nonlinear operator variable');
end

fcn = extstruct.fcn;
try
    n = yalmip('nvars');
    [F,properties,arguments] = feval(fcn,method,extstruct.var,extstruct.arg{1:end-1});
    if isa(F,'constraint')
        F = lmi(F);
    end
    newAux = n+1:yalmip('nvars');
    involved = getvariables(extstruct.arg{1});
    for i = 2:length(extstruct.arg)-1
        vars = getvariables(extstruct.arg{i});
        if ~isempty(vars)           
            involved = union(involved,vars);
        end
    end
    if ~isempty(options)
        if ~(strcmp(options.robust.auxreduce,'none'))
            % This info is only needed when we do advanced Robust optimization
            yalmip('setdependence',[getvariables(extstruct.var) newAux],involved);
            yalmip('setdependence',[getvariables(extstruct.var)],newAux);
        end
    end
catch
    error(['Failed when trying to create a model for the "' extstruct.fcn '" operator']);
end

% Make sure all operators have these properties, and that they are
% consistently set (positive has non-negative domain etc)
if ~isempty(properties)
    if ~iscell(properties)
        properties = {properties};
    end
    for i = 1:length(properties)
        properties{i}.name = fcn;
        properties{i} = assertOperatorProperties(properties{i});    
        properties{i} = assertProperty(properties{i},'models',getvariables(extstruct.var));
    end
end

% Normalize the callback expression if necessary, add domain constraints, 
% and check for some obsoleted stuff
if ~isempty(properties)
    if isequal(properties{1}.model,'callback')
        [F_normalizing,normalizer,argument] = NormalizeCallback(method,extstruct.var,extstruct.arg{:},options.warmstart);
        F = F + F_normalizing;
        if ~isempty(normalizer)
            if ~isempty(properties{1}.domain)
                if ~isinf(properties{1}.domain(1))
                    F = [F, normalizer >= properties{1}.domain(1)];
                end  
                if ~isinf(properties{1}.domain(2))
                    F = [F, normalizer <= properties{1}.domain(2)];
                end
            end
        end
        if ~isempty(argument)
            if ~isempty(properties{1}.domain)
                if ~isinf(properties{1}.domain(1))
                    F = [F, argument >= properties{1}.domain(1)];
                end  
                if ~isinf(properties{1}.domain(2))
                    F = [F, argument <= properties{1}.domain(2)];
                end
            end
        end
    end
    if length(extstruct.computes)>1
        for i = 1:length(properties)
            properties{i}.models = extstruct.computes;
        end
    end
    for i = 1:length(properties)
        if ~isa(properties{i}.convexity,'function_handle')
            if ~any(strcmpi(properties{i}.convexity,{'convex','concave','none'}))
                disp('More cleaning, strange convextiy returned...Report bug in model.m')
                error('More cleaning, strange convextiy returned...Report bug in model.m')
            end
        end
    end
end

% This is useful in MPT
if ~isempty(F)
    F = tag(F,['Expansion of ' extstruct.fcn]);
end
