function F = le(X,Y)
% Internal class for constraint lists

% Try to evaluate
try
    if isa(X,'constraint')
        % (z > w) < y
        Z = Y - X.List{end};
        F = X;
        F.List{end+1} = '<=';
        F.List{end+1} = Y;
        F.Evaluated{end+1} = Z;
        F.ConstraintID(end+1) = yalmip('ConstraintID');
        F.strict(end+1) = 0;
    else
        % x < (w > y)
        Z = Y.List{1} - X;
        F = Y;
        F.List = {X,'<=',F.List{:}};
        F.Evaluated = {Z,F.Evaluated{:}};
        F.ConstraintID = [yalmip('ConstraintID') F.ConstraintID];
        F.strict = [1 F.strict];
    end
catch
    error(lasterr);
end


