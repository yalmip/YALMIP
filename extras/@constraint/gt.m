function F = gt(X,Y)
% Internal class for constraint lists

superiorto('sdpvar');
superiorto('double');

warned = 0;
if isa(X,'constraint')    
    if ~is(sdpvar(X),'quantized')    
       warning('YALMIP:strict','Strict inequalities will only (if at all) be supported on binvar/intvar variables in future versions of YALMIP. Turn off this warning using warning(''off'',''YALMIP:strict'')');
       %error('Non-strict inequalities will only be supported on binvar/intvar variables in future versions of YALMIP');
    end
end
if warned == 0
    if isa(Y,'constraint')        
        if ~is(sdpvar(Y),'quantized')        
           warning('YALMIP:strict','Strict inequalities will only (if at all) be supported on binvar/intvar variables in future versions of YALMIP. Turn off this warning using warning(''off'',''YALMIP:strict'')');
           %error('Non-strict inequalities will only be supported on binvar/intvar variables in future versions of YALMIP');
        end
    end
end

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
        F.strict(end+1) = 1;
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

