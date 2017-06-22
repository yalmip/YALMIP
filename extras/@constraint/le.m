function F = le(X,Y)
% Internal class for constraint lists

% Try to evaluate
try
    if isa(X,'constraint')
        % (z > w) < y
        try        
            Z = Y - X.List{end};
        catch
            Y = reshape(Y,[],1);
            Z = Y - X.List{end};
        end
        F = X;
        F.List{end+1} = '<=';
        F.List{end+1} = Y;
        F.Evaluated{end+1} = Z;
        F.ConstraintID(end+1) = yalmip('ConstraintID');
        F.strict(end+1) = 0;
    else
        % x < (w > y)
        try
            Z = Y.List{1} - X;
        catch
            X = reshape(X,[],1);
            Z = Y.List{1} - X;
        end
        F = Y;
        F.List = {X,'<=',F.List{:}};
        F.Evaluated = {Z,F.Evaluated{:}};
        F.ConstraintID = [yalmip('ConstraintID') F.ConstraintID];
        F.strict = [1 F.strict];
    end
catch
    error(lasterr);
end


