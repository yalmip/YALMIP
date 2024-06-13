function y = eq(X,Y)
%EQ (overloaded)

if isa(X,'blkvar')
    X = sdpvar(X);
end

if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

try
    Z = X - Y;
    if ishermitian(Z)
    	if isreal(Z)
		    mask = triu(true(length(Z)));
		    y = constraint(extsubsref(Z,mask),'==', 0);
		else
		    mask = triu(true(length(Z)),1);
    	    y = constraint([diag(Z);real(extsubsref(Z,mask));imag(extsubsref(Z,mask))],'==', 0);
    	end
    else
        y = constraint(Z,'==',0);
    end
catch
    error(lasterr)
end
