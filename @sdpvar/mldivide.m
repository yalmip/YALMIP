function y = mldivide(X,Y)
%MLDIVIDE (overloaded)

if (isa(X,'sdpvar'))
  error('Division of matrix variables not possible.')
end

try
  lmi_variables = getvariables(Y);
  nv = length(lmi_variables);
  y  = Y;
  n = Y.dim(1);
  m = Y.dim(2);
  if m==1     
    y.basis = X\Y.basis;   
    y.dim(1) = size(y.basis,1);
    y.dim(2) = 1;
  else % FIX : VECTORIZE THIS...
    [L,U] = qr(X);
    temp = U\(L\reshape(Y.basis(:,1),n,m));
    y.basis = temp(:);
    for i = 1:nv
        temp = U\(L\reshape(Y.basis(:,i+1),n,m));
        y.basis(:,i+1) = temp(:);
    end;
    y.dim(1) = size(temp,1);
    y.dim(2) = size(temp,2);
  end
  y = flush(y);
catch
  error(lasterr);
end