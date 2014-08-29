function Q=getbasevectorwithoutcheck(X,ind)
%GETBASEVECTORWITHOUTCHECK Internal function to extract basematrix for variable ind

Q=X.basis(:,ind+1);
  
  
      