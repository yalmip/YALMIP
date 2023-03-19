function y = eq(X,Y)
%EQ (overloaded)

if isa(X,'blkvar')
    X = sdpvar(X);
end

if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

try
    if ishermitian(X) && ishermitian(Y)
    	if isreal(X) && isreal(Y)
		    mask = triu(true(length(X)));
		    y = constraint(extsubsref(X,mask),'==',extsubsref(Y,mask));
		else
		    mask = triu(true(length(X)),1);
    	    y = constraint([diag(X);real(extsubsref(X,mask));imag(extsubsref(X,mask))],'==',[diag(Y);real(extsubsref(Y,mask));imag(extsubsref(Y,mask))]);
    	end
    else
        y = constraint(X,'==',Y);
    end
catch
    error(lasterr)
end
