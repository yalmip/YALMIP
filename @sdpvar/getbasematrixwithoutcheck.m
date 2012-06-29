function Q=getbasematrixwithoutcheck(X,ind)
%GETBASEMATRIXWITHOUTCHECK Internal function to extract basematrix for variable IND

% Author Johan Löfberg 
% $Id: getbasematrixwithoutcheck.m,v 1.3 2006-07-26 20:17:58 joloef Exp $  

Q=reshape(X.basis(:,ind+1),X.dim(1),X.dim(2));
  
  
      