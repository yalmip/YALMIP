function Q=getbase(X,z)
%GETBASE Internal function to extract all base matrices

if nargin == 1
    Q=X.basis;
else
    Q = [];
    for i = 1:length(z)
        q = getbasematrix(X,z.lmi_variables(i));
        Q = [Q q(:)];
    end
end
  
  
      