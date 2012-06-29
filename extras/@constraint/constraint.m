function F = constraint(X,quantifier,Y)
% Internal class for constraint list

% Author Johan Löfberg
% $Id: constraint.m,v 1.8 2009-04-29 07:48:12 joloef Exp $

superiorto('sdpvar');
superiorto('double');
superiorto('logical');

if isa(X,'blkvar')
    X = sdpvar(X);
end
if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

% Evaluate
switch quantifier
    case '>'
        Z = X - Y;      
    case '>='
        Z = X - Y;     
    case '<'
        Z = Y - X;      
    case '<='
        Z = Y - X;       
    case '=='
        Z = Y - X;      
    otherwise
        error('Quantifier not supported')
end

if isa(Z,'double')
    
    if  size(Z,1)==size(Z,2) &&  norm(Z-Z',inf)<1e-12 && ~isequal(quantifier,'==')
        checkSDP = 1;
    else
        checkSDP = 0;
    end
    
    if checkSDP
        if min(eig(Z))>=0
            warning('Constraint evaluated to trivial true.')
            F = [];
            return
        else
            error('SDP constraint evaluated to trivial false (no decision variable in constraint)')            
        end
    else
        Z = Z(:);
        switch quantifier
            case '=='
                if all(Z)==0
                    warning('Constraint evaluated to trivial true.')
                    F = [];
                    return
                else
                    error('Equality constraint evaluated to trivial false (no decision variable in constraint)')
                end
            case {'<=','>='}
                if all(Z>=0)
                    warning('Constraint evaluated to trivial true.')
                    F = [];
                    return
                else
                    error('Inequality constraint evaluated to trivial false (no decision variable in constraint)')
                end
            case {'<','>'}
                if all(Z>0)
                    warning('Constraint evaluated to trivial true.')
                    F = [];
                    return
                else
                    error('Inequality constraint evaluated to trivial false (no decision variable in constraint)')
                end
        end
    end
end

switch quantifier
case {'>','<'}
    F.strict(1) = 1;
case {'>=','<=','=='}
    F.strict(1) = 0;
otherwise
    error('Quantifier not supported')
end

F.List={X,quantifier,Y};
F.Evaluated{1} = Z;
F.ConstraintID = yalmip('ConstraintID');
F.tag{1} = '';
F = class(F,'constraint');
	