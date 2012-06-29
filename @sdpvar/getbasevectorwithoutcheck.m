function Q=getbasevectorwithoutcheck(X,ind)
%GETBASEVECTORWITHOUTCHECK Internal function to extract basematrix for variable ind

% Author Johan Löfberg 
% $Id: getbasevectorwithoutcheck.m,v 1.2 2004-07-01 11:17:10 johanl Exp $  

Q=X.basis(:,ind+1);
  
  
      