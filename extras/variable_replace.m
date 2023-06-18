function Z = variable_replace(X,Y,W)

[mt,variabletype] = yalmip('monomtable');
if isequal(getbase(Y),[0 1]) && isequal(getbase(W),[0 1])    
    % Super standard case, replace monom with other monom
    % However, it is a nonlinear monom to be replaced
    % so it must mean monom replacement, not actually replacing a variable
    % i.e. replace(x+x^2,x,y) means y+y^2 but we have
    % i.e. replace(x+x^2,x^2,y) meaning x+y but we have
    % Special case added for the nonlinear kkt rewrite to work
    k = getvariables(Y);
    if variabletype(k)
       % Yep, that's the case here
       k_new = getvariables(W);
       Z = X;
       Z = struct(Z);
       i = find(Z.lmi_variables == k);
       if ~isempty(i)         
        Z.lmi_variables(i) = k_new;
        % Messed up order (lmi_variables should be sorted)
        if any(diff(Z.lmi_variables)<0)
            [i,j]=sort(Z.lmi_variables);
            Z.basis = [Z.basis(:,1) Z.basis(:,j+1)];
            Z.lmi_variables = Z.lmi_variables(j);
        end
        Z = sdpvar(Z.dim(1),Z.dim(2),[],Z.lmi_variables,Z.basis);
       else
           Z = X;
       end          
       return
    end
end

% Check so that Y is a simple unit variable
Ybase = getbase(Y);
Yvariables = getvariablesSORTED(Y);
Xbase = getbase(X);
Xvariables = getvariables(X);
[i,j,k] = find(Ybase);
if ~isequal(sort(i),1:length(i))
end
if ~isequal(sort(j),2:(length(i)+1))
end
if ~all(k == 1)
end

% Linear, or at least linear in Y
if all(variabletype(Xvariables) == 0) %| all(sum(mt(any(mt(getvariables(X),getvariables(Y)),2),:),2)==1)
    % Simple linear replacement    
    [aa,bb] = ismember(Xvariables,Yvariables);
    i1 = find(aa);
    v1 = bb(i1);
    i2 = find(~aa);
    v2 = Xvariables(i2);    
    v = sparse(i1,ones(length(i1),1),W(v1),length(Xvariables),1);
    v = v + sparse(i2,ones(length(i2),1),recover(v2),length(Xvariables),1);    
    Z = Xbase*[1;v];    
    Z = reshape(Z,size(X,1),size(X,2));
else
    if nnz(mt(getvariables(X),getvariables(Y)))==0
        Z = X;
    else
        Z = nonlinearreplace(X,Y,W);
    end
    return
end

% This has not been tested (copied from variable_replace) so it is placed
% in a catch to be safe.
try
    Xvariables = getvariables(Z);
    extvar = yalmip('extvariables');
    Xext = find(ismember(Xvariables,extvar));
    if ~isempty(Xext)
        %We must dig down in extended operators to see if they use the replaced
        %set of variables
        for i = 1:length(Xext)
            extstruct = yalmip('extstruct',Xvariables(Xext(i)));
            anychanged = 0;
            for j = 1:length(extstruct.arg)
                if isa(extstruct.arg{j},'sdpvar')
                    XinY = find(ismembc(getvariables(extstruct.arg{j}),Yvariables));
                    if ~isempty(XinY)
                        anychanged = 1;
                        extstruct.arg{j} = replace(extstruct.arg{j},Y,W);
                    else
                    end
                end
            end
            if anychanged
                Zi = yalmip('define',extstruct.fcn,extstruct.arg{:});
                Xvariables(Xext(i)) = getvariables(Zi);
            end
        end
        % And now recover this sucker
        Z = struct(Z);
        Z.lmi_variables = Xvariables;
        % Messed up order (lmi_variables should be sorted)
        if any(diff(Z.lmi_variables)<0)
            [i,j]=sort(Z.lmi_variables);
            Z.basis = [Z.basis(:,1) Z.basis(:,j+1)];
            Z.lmi_variables = Z.lmi_variables(j);
        end
        Z = sdpvar(Z.dim(1),Z.dim(2),[],Z.lmi_variables,Z.basis);
    end
catch
end


function Yvariables = getvariablesSORTED(Y)
B = getbase(Y);
if any(B(:,1))
  error
end
B = B(:,2:end);
if ~(all(sum(B,2)==1) && all(sum(B,2)==1))
 	error
 end
v = getvariables(Y);
[i,j,k] = find(B);
Yvariables = v(i);

