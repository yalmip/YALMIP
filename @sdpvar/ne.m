function F = ne(X,Y)
%NE (overloaded)
%
%    F = ne(x,y)
%
%   See also SDPVAR/AND, SDPVAR/OR, BINVAR, BINARY

% Models NE using logic constraints

% bin1 = isa(X,'sdpvar') | isa(X,'double');
% bin2 = isa(Y,'sdpvar') | isa(Y,'double');
%
% if ~(bin1 & bin2)
%     error('Not equal can only be applied to integer data')
% end

if isa(X,'sdpvar') & isa(Y,'sdpvar') & is(X,'binary') & is(Y,'binary')
    F = (X + Y == 1);
elseif isa(X,'sdpvar') & is(X,'binary') & is(X,'lpcone') &  ~isa(Y,'sdpvar') &  ismember(Y,[0 1])
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
elseif isa(X,'sdpvar') & is(X,'integer')% &  isa(Y,'double')
    X = reshape(X,[],1);
    Y = reshape(Y,[],1);
    if length(X)>1 & length(Y)==1
        Y = repmat(Y,length(X),1);
    end
    less = binvar(length(X),1);
    larger = binvar(length(X),1);
    F = [implies(less,X<=Y-1), implies(larger,X>=Y+1), sum(less) + sum(larger) >= 1];
else
    F = ((X<=Y-0.5) | (X>=Y+0.5));
end