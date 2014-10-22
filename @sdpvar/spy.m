function varargout = spy(X)
%SPY (overloaded)

 if isa(X,'blkvar')
    X = sdpvar(X);
 end
    
Z = reshape(sum(abs(X.basis),2),X.dim(1),X.dim(2));
Z = Z~=0;
if nargout==0
    spy(Z)
else
    varargout{1}=Z;
end