function varargout = ismember(varargin)
%ISMEMBER Define membership constraint on SDPVAR object
%
% F = ISMEMBER(x,P)
%
% Input
%  x : SDPVAR object 
%  P : MPT polytope object, double, or CONSTRAINT object
% Output
%  F : Constraints
%
% Depending on the second argument P, different classes of constraint are
% generated. 
%
% If P is a single polytope, the linear constraints [H,K] = double(P);
% F=[H*x <= K] will be created. 
%
% If P is a polytope array, then length(P) binary variables will be
% introduced and the constraint will model that x is inside at least one of
% the polytopes. 
%
% If P is a vector of DOUBLE, a constraint constraining the elements of x
% to take one of the values in P is created. This will introduce
% numel(P)*numel(x) binary variables 
%
% If P is matrix of DOUBLE, a constraint constraining the vector x to equal
% one of the columns of P is created. This will introduce size(P,2) binary
% variables  
%
% Since the two last constructions are based on big-M formulations, all
% involved variable should have explicit variable bounds. 

x = varargin{1};
p = varargin{2};

% Backwards compatibility (this should really be done in another command)
% This code is probably only used in solvemoment
if isa(x,'double')
    varargout{1} = any(full(p.basis(:,1)));
    return
end

if isa(x,'sdpvar') & isa(p,'sdpvar')
    if numel(x) > 1
        d = size(x);
        x = reshape(x,[],1);
        YESNO = [];
        for i = 1:length(x)
            YESNO = [YESNO;ismember(subsref(x,struct('type','()','subs',{{i}})),p)];
        end
        varargout{1} = reshape(YESNO,d);        
    else
        d = size(p);
        p = reshape(p,[],1);
        varargout{1} = false;
        for k = 1:length(p)
            pk = subsref(p,struct('type','()','subs',{{k}}));
            if isequal(getbase(pk),getbase(x)) && isequal(getvariables(pk),getvariables(x))
                varargout{1} = true;
                return 
            end
        end        
    end
    return
end

% Here is the real overloaded ismember
switch class(varargin{1})
    case 'sdpvar'
        if isa(varargin{1},'sdpvar') & (isa(varargin{2},'polytope') | isa(varargin{2},'Polyhedron'))
            if ~isequal(length(varargin{1}),safe_dimension(varargin{2}))
                disp('The polytope in the ismember condition has wrong dimension')
                error('Dimension mismatch.');
            end
        end
        if isa(varargin{2},'polytope') & length(varargin{2})==1
            [H,K] = double(varargin{2});
            varargout{1} = [H*x <= K];
        elseif isa(varargin{2},'Polyhedron') & length(varargin{2})==1
            %P = convexHull(varargin{2});
            P = minHRep(varargin{2});
            varargout{1} = [P.A*x <= P.b, P.Ae*x == P.be];
        else                            
            varargout{1} = (yalmip('define',mfilename,varargin{:}) == 1);            
            varargout{1} = setupMeta(lmi([]), mfilename,varargin{:});
            if isa(varargin{2},'double')
                if size(varargin{1},1) == size(varargin{2},1) % v in [v1 v2 v3]
                 varargout{1} = [ varargout{1}, min(varargin{2},[],2) <= varargin{1} <= max(varargin{2},[],2)];
                else
                 varargout{1} = [ varargout{1}, min(min(varargin{2})) <= varargin{1}(:) <= max(max(varargin{2}))];
                end
            end
        end

    case 'char'
        varargout{1} = ismember_internal(varargin{3},varargin{4});
end

function d = safe_dimension(P)
if isa(P,'polytope')
    d = dimension(P);
elseif isa(P,'Polyhedron')
    d = P.Dim;
end