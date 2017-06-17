function F = exclude(X,Y)
%EXCLUDE Excludes a binary solution
%
%    F = exclude(X,value)
%
%EXCLUDE is used to avoid a particular binary solution. This can be used
% to repeatedly solve MILP problems while exluding all past solutions
%
% A = randn(30,15);
% b = 25*rand(30,1);
% c = randn(15,1);
% x = binvar(15,1);
% Model = A*x <= b;
% sol = solvesdp(Model,c'*x);
% while sol.problem == 0
%    Model = [Model, exclude(x,double(x))];
%    sol = solvesdp(Model,c'*x);
% end

if isa(X,'sdpvar') & is(X,'binary') &  isnumeric(Y) &  ismember(Y,[0 1])
    
    if isequal(size(X),size(Y))
    else
        error('Dimension mismatch in EXCLUDE')
    end
    
    zv = find((Y == 0));
    ov = find((Y == 1));
    lhs = 0;
    if ~isempty(zv)
        lhs = lhs + sum(extsubsref(X,zv));
    end
    if ~isempty(ov)
        lhs = lhs + sum(1-extsubsref(X,ov));
    end
    F = [lhs >=1];
    
else
    error('EXCLUDE only applicable to binary variables and data');
end