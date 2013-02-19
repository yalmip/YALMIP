function F = exclude(X,Y)
%EXCLUDE Excludes a binary solution
%
%    F = exclude(X,value)
%
% EXCLUDE is used to avoid a particular binary solution. This can be used
% to repeatedly solve MILP problems while exluding all past solutions.

if isa(X,'sdpvar') & is(X,'binary') &  isa(Y,'double') &  ismember(Y,[0 1])
    
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