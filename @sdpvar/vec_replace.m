function Z = replace(X,Y,W,expand)
%REPLACE Substitutes variables
%
%Z = REPLACE(Y,X,W)  Replaces any occurence of the SDPVAR object Y
%                    in the SDPVAR object X with the double W
%
% Example
%  x = sdpvar(1,1);
%  t = sdpvar(1,1);
%  Y = [1+t;1+x+t];
%  Y = replace(Y,x,2) generates Y=[1+t;3+t]

% if nargin<4
%     expand = 1;
% end
% 
% if ~isa(X,'sdpvar')
%     Z = X;
%     return
% end
% if ~isa(Y,'sdpvar')
%     error('Second arguments must be an sdpvar object')
% end
% 
% if ~is(Y,'linear') 
%     error('Second arguments must be linear')
% end
% 
% if prod(size(W)) == 1
%     W = repmat(W,size(Y));
% end
% 
% if ~isequal(size(Y),size(W))
%     if isequal(fliplr(size(Y)),size(W))
%         W = W';
%     else
%     error('Both arguments must have same size')
%     end
% end

% if isa(W,'sdpvar')
%     % This is tricky...
%     Z = variable_replace(X,Y,W);
%     return
% end
% 
% if ~isa(W,'double')
%     error('Third arguments must be a double')
% end

% Replace with NaN   destroys everything, assume it should be cleared
W(isnan(W)) = 0;

y_lmi_variables = Y.lmi_variables;
% b = W(:)-Y.basis(:,1);
% A = Y.basis(:,2:end);
% feas_var = A\b;
% if norm(A*feas_var-b)>sqrt(eps)
%     error('Inconsistent assignment')
% end

x_lmi_variables = X.lmi_variables;
n = X.dim(1);
m = X.dim(2);

[monomtable,variabletype] = yalmip('monomtable');
if all(variabletype(x_lmi_variables)==0) % is(X,'linear')
    Z = X.basis(:,1);
    for i = 1:length(x_lmi_variables)
        j = find(x_lmi_variables(i) == y_lmi_variables);
        if isempty(j)
            Z = Z + recover(x_lmi_variables(i))*X.basis(:,i+1);
        else
            Z = Z + feas_var(j)*X.basis(:,i+1);
        end
    end
else
    replaced_vars = depends(Y);
    % used_variables = getvariables(X);
    used_variables = x_lmi_variables;
    %  monomtable = yalmip('monomtable');
    local_monom = monomtable(used_variables,replaced_vars);
    Wall = W;
    Z = [];

    local_monoms_left = monomtable(used_variables,:);
    local_monoms_left(:,replaced_vars) = 0;
    used_left = find(sum(local_monoms_left,1));
    base = recovermonoms(local_monoms_left(:,used_left),recover(used_left));
                
    for rep = 1:size(Wall,2)
        W = Wall(:,rep);
        W = W(:)';
        gain = zeros(length(used_variables),1);
        for i = 1:length(used_variables)
            % F**N 6.5 0^sparse(0) and 0^0 differ
            gain(i) = prod(W.^full(local_monom(i,:)));
        end

        base = base.*gain(:);
        Ztemp = X.basis(:,1);
        Ztemp = Ztemp + X.basis(:,2:end)*base;
        Z = [Z;Ztemp];
    end
end

if isa(Z,'sdpvar')
    Z.dim(1) = n*size(Wall,2);
    Z.dim(2) = m;
    Z.typeflag = X.typeflag;
    % Reset info about conic terms
    Z.conicinfo = [0 0];
else
    Z = reshape(full(Z),n,m);
end
