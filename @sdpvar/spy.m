function varargout = spy(X)
%SPY (overloaded)

% Author Johan Löfberg 
% $Id: spy.m,v 1.5 2006-07-26 20:17:58 joloef Exp $   

 if isa(X,'blkvar')
    X = sdpvar(X);
 end
    
Z = reshape(sum(abs(X.basis),2),X.dim(1),X.dim(2));
if nargout==0
    spy(Z)
else
    varargout{1}=Z;
end