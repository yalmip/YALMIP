function y = gt(X,Y)
%GT (overloaded)

% Author Johan Löfberg
% $Id: gt.m,v 1.3 2005-06-17 13:02:01 joloef Exp $

if isa(X,'blkvar')
    X = sdpvar(X);
end

if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

warned = 0;
if isa(X,'sdpvar')
    if ~is(X,'quantized')    
       warning('YALMIP:strict','Strict inequalities will only (if at all) be supported on binvar/intvar variables in future versions of YALMIP. Turn off this warning using warning(''off'',''YALMIP:strict'')');
       %error('Non-strict inequalities will only be supported on binvar/intvar variables in future versions of YALMIP');
    end
end
if warned == 0
    if isa(Y,'sdpvar')
        if ~is(Y,'quantized')        
           warning('YALMIP:strict','Strict inequalities will only (if at all) be supported on binvar/intvar variables in future versions of YALMIP. Turn off this warning using warning(''off'',''YALMIP:strict'')');
           %error('Non-strict inequalities will only be supported on binvar/intvar variables in future versions of YALMIP');
        end
    end
end

try
    y = constraint(X,'>',Y);
catch
    error(lasterr)
end
