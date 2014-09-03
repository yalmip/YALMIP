function [F,properties,arguments,fcn]=model(X,method,options,extstruct)
%MODEL  Extracts nonlinear operator models
%
% [F,properties] = model(x)
%
% MODEL returns the constraints needed to model a variable related to an
% extended operator such as min, max, abs, norm, geomean, ...
%
% Examples :
%
% sdpvar x y;
% t = min(x,y);
% [F,properties] = epigraph(t)
% Gives (F = (t<=x) + (t<=y))
%
% sdpvar x y
% t = max(norm([x;y],1+y))
% [F,properties] = epigraph(t)
% Gives (F = (u<=t) + (1+y<=t))
% where u is the variable modelling norm([x;y])

extvar = getvariables(X);
arguments   = cell(1,length(extvar));
properties  = cell(1,length(extvar));

if nargin<2
    method = 'graph';
end

if nargin < 3
    options = [];
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
switch fcn
    % *********************************************************************
    % mpower, max, min  are currently implmeneted in a non-standard way for
    % both performance and technical reasons
    % *********************************************************************
    case 'mpower'
        [F,properties,arguments] = mpower_internal(X,method,options,extstruct);

    case 'max'
        [F,properties,arguments] = max_internal(X,method,options,extstruct);

    case 'min'
        [F,properties,arguments] = min_internal(X,method,options,extstruct);

    otherwise
        try
            [F,properties,arguments] = feval(extstruct.fcn,method,extstruct.var,extstruct.arg{:});
        catch
            error(['Failed when trying to create a model for the "' extstruct.fcn '" operator']);
        end
end
% This field is not official, and currently only used in sort
if ~isfield(properties,'models')
    properties.models = getvariables(extstruct.var);
end


