function A=accumarray(subs,val,varargin)
%ACCUMARRAY (overloaded)

if size(subs,2)>1 || min(size(val))>1 || nargin > 2
    error('SDPVAR/ACCUMARRAY currently only supports simple cases (subs=column, val=vector)');
end

S = sparse(max(subs),length(val),0);
for i = 1:size(S,1)
    S(i,find(i == subs)) = 1;
end
A = S*reshape(val,[],1);