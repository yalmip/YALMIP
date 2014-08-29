function y = kron(X,Y)
%KRON (overloaded)

if isequal(X,1)
    y = Y;
    return;
end

if isequal(Y,1);
    y = X;
end

if (isa(X,'sdpvar') & isa(Y,'sdpvar'))
    y = [];
    for i = 1:size(X,1)
        this_row = [];
        s.type = '()';
        for j = 1:size(X,2);
            s.subs = {i,j};
            this_row = [this_row subsref(X,s)*Y];
        end
        y = [y;this_row];
    end
    return
end

if isa(X,'sdpvar')       
	lmi_variables = getvariables(X);
	nv = length(lmi_variables);
	y  = X;
	sparse_Y = sparse(Y);
	temp = kron(getbasematrix(X,0),sparse_Y);
	y.basis=temp(:);
	for i = 1:nv
		temp = kron(getbasematrix(X,lmi_variables(i)),sparse_Y);
		y.basis(:,i+1) = temp(:);
	end;	
end

if isa(Y,'sdpvar')
  lmi_variables = getvariables(Y);
  nv = length(lmi_variables);
  y  = Y;
  sparse_X = sparse(X);
  temp = kron(sparse_X,getbasematrix(Y,0));  
  y.basis = temp(:); 
  for i = 1:nv
		temp = kron(sparse_X,getbasematrix(Y,lmi_variables(i)));
    y.basis(:,i+1) = temp(:);
  end;
end
y.dim(1) = size(temp,1);
y.dim(2) = size(temp,2);
y = clean(y);
% Reset info about conic terms
y.conicinfo = [0 0];