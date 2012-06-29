function F = gt(X,Y)
% Internal class for constraint lists

% Author Johan Löfberg
% $Id: ge.m,v 1.2 2007-09-12 14:28:29 joloef Exp $

superiorto('sdpvar');
superiorto('double');

try
    % Try to evaluate
    if isa(X,'constraint')
        % (z > w) > y
        Z = X.List{end} - Y;
        F = X;
        F.List{end+1} = '>';
        F.List{end+1} = Y;
        F.Evaluated{end+1} = Z;
        F.ConstraintID(end+1) = yalmip('ConstraintID');
        F.strict(end+1) = 0;
    else
        % z > (w > y)
        Z = X - Y.List{1};
        F = Y;
        F.List = {X,'>',F.List{:}};
        F.Evaluated = {Z,F.Evaluated{:}};
        F.ConstraintID = [yalmip('ConstraintID') F.ConstraintID];        
        F.strict = [1 F.strict];               
    end
catch
    error(lasterr);
end

