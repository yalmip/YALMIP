function varargout = getvariables(X,whatvars)
%GETVARIABLES Returns variable indicies to variables used in a SDPVAR object
%
% i = getvariables(x)
%
% Input
%    x : SDPVAR object
% Output
%    i : DOUBLE
%
% Note : The difference compared DEPENDS is that nonlinear terms
%        are reported as separate variables.

% Author Johan Löfberg
% $Id: getvariables.m,v 1.3 2004-10-04 07:51:13 johanl Exp $

switch nargin
    case 1
        varargout{1} = X.lmi_variables;
    case 2
        [sqrList,CompressedList] = yalmip('nonlinearvariables');
        used_variables = getvariables(X);
        if isempty(sqrList)
            LinearVariables = used_variables;
            NonlinearVariables=[];
        else
            LinearVariables    = setdiff(used_variables,sqrList(:,1));
            NonlinearVariables = setdiff(used_variables,LinearVariables);
        end
        switch whatvars
            case 'linear'
                varargout{1} = LinearVariables;
            case 'nonlinear'
                varargout{1} = NonlinearVariables;
            case 'both'
                varargout{1} = LinearVariables;
                varargout{2} = NonlinearVariables;
            otherwise
                error
        end
    otherwise
        error
end
