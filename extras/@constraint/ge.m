function F = ge(X,Y)
% Internal class for constraint lists

try
    % Try to evaluate
    if isa(X,'constraint')
        % (z > w) > y
        try
            Z = X.List{end} - Y;
        catch
            Y = reshape(Y,[],1);
            Z = X.List{end} - Y;
        end
        F = X;
        F.List{end+1} = '>=';
        F.List{end+1} = Y;
        F.Evaluated{end+1} = Z;
        F.ConstraintID(end+1) = yalmip('ConstraintID');
        F.strict(end+1) = 0;
    else
        % z > (w > y)
        try
            Z = X - Y.List{1};
        catch
            X = reshape(X,[],1);
            Z = X - Y.List{1};
        end
        F = Y;
        F.List = {X,'>=',F.List{:}};
        F.Evaluated = {Z,F.Evaluated{:}};
        F.ConstraintID = [yalmip('ConstraintID') F.ConstraintID];        
        F.strict = [1 F.strict];               
    end
catch
    error(lasterr);
end

