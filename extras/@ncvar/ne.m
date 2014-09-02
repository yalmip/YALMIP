function F = ne(X,Y)
%NE (overloaded)
%
%    F = (ne(x,y))
%
%   See also SDPVAR/AND, SDPVAR/OR, BINVAR, BINARY

% Models NE using logic constraints

% bin1 = isa(X,'sdpvar') | isa(X,'double');
% bin2 = isa(Y,'sdpvar') | isa(Y,'double');
%
% if ~(bin1 & bin2)
%     error('Not equal can only be applied to integer data')
% end

if is(X,'binary') &  isa(Y,'double') & all((Y == round(Y)))
    zv = find((Y == 0));
    ov = find((Y == 1));
    lhs = 0;
    if ~isempty(zv)
        lhs = lhs + sum(extsubsref(X,zv));
    end
    if ~isempty(ov)
        lhs = lhs + sum(1-extsubsref(X,ov));
    end
    F = (lhs >=1);
else
    F = ((X<=Y-0.5) | (X>=Y+0.5));
end