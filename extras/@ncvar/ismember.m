function varargout = ismember(varargin)
%ISMEMBER Define membership constraint on SDPVAR object
%
% F = ISMEMBER(x,P)
%
% Input
%  x : SDPVAR object 
%  P : MPT polytope object or double
% Output
%  F : SET object
%
% Depending on the second argument P, different classes of constraint are
% generated. 
%
% If P is a single polytope, the linear constraints [H,K] = double(P);
% F=set(H*x <= K) will be created. 
%
% If P is a polytope array, then length(P) binary variables will be
% introduced and the constraint will model that x is inside at least one of
% the polytopes. 
%
% If P is a DOUBLE, a constraint constraining the elements of x to take one
% of the values in P is created. This will introduce numel(P)*numel(x)
% binary variables 
%
% Since the two last constructions are based on big-M formulations, all
% involved variable should have explicit variable bounds. 

% Author Johan Löfberg 
% $Id: ismember.m,v 1.1 2006-08-10 18:00:20 joloef Exp $  


x = varargin{1};
p = varargin{2};

% Backwards compatibility (this should really be done in another command)
% This code is probably only used in solvemoment
if isa(x,'double')
    varargout{1} = any(full(p.basis(:,1)));
    return
end

if isa(x,'sdpvar') & isa(p,'sdpvar')
    x_base = x.basis;
    x_vars = x.lmi_variables;

    p_base = x.basis;
    p_vars = x.lmi_variables;

    % Member at all
    varargout{1} = ismember(x.lmi_variables,p.lmi_variables);
    if varargout{1}
        index_in_x_vars = find(x.lmi_variables == p.lmi_variables);
        varargout{1} = full(any(p.basis(:,1+index_in_x_vars),2));
        if min(p.dim(1),p.dim(2))~=1
            varargout{1} = reshape(YESNO,p.dim(1),p.dim(2));
        end
    end
    return
end

% Here is the real overloaded ismember
switch class(varargin{1})
    case 'sdpvar'
        varargout{1} = set(yalmip('addextendedvariable',mfilename,varargin{:}) == 1);

    case 'char'
        varargout{1} = ismember_internal(varargin{3},varargin{4});
        varargout{2} = struct('convexity','milp','monotonicity','milp','definiteness','positive','extra','marker');
        varargout{3} = varargin{3};
end
