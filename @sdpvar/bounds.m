function varargout = bounds(x,lower,upper,aux)
%BOUNDS Adds implicit bounds on variables.
%
% BOUNDS IS OBSOLETE: Use standard SET constraints
%
% BOUNDS(x,Lower,Upper)   Adds bound constraints on variables.
%                         These bounds are used when performing
%                         big M formulations and similar things.

% Author Johan Löfberg 
% $Id: bounds.m,v 1.6 2006-05-11 10:49:13 joloef Exp $   

variables = getvariables(x);
if nargin == 1
    lb = yalmip('getbounds',variables);
    lower = lb(:,1);
    upper = lb(:,2);
else
    lower = lower(:);
    upper = upper(:);
    if length(lower)==1
        lower = repmat(lower,length(variables),1);
    end
    if length(upper)==1
        upper = repmat(upper,length(variables),1);
    end
    % 0 - No problems
    % 1 - Trying to bound nonlinear variable
    % 2 - Trying to bound nonlinear operator     
    fail = yalmip('setbounds',variables,lower,upper);
    switch fail
        case {1,2}
            error('BOUNDS can only be applied to linear unitary SDPVAR variables.')
        otherwise
    end
end

if nargout>0
    varargout{1} = lower;
    varargout{2} = upper;
end