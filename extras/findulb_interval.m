function [lb,ub,cand_rows] = findulb(F_struc,K)
%FINDULB Internal function to extract upper and lower variable bounds

% Author Johan Löfberg
% $Id: findulb_interval.m,v 1.1 2006-12-14 15:37:52 joloef Exp $

% special code for the interval data case, to avoid overhead in the default
% function findulb 

n = size(F_struc,2)-1;
lb = -inf*ones(n,1);
ub = inf*ones(n,1);
cand_rows = [];

if (K.f ~=0)
    A = -F_struc(1:K.f,2:end);
    b = F_struc(1:K.f,1);
    n = size(F_struc,2)-1;
    cand_rows = find(sum(A~=0,2)==1);
    for i = 1:length(cand_rows)
        j = find(A(cand_rows(i),:));
        ub(j)=min(ub(j),b(cand_rows(i))/A(cand_rows(i),j));
        lb(j)=max(lb(j),b(cand_rows(i))/A(cand_rows(i),j));
    end
end

if (K.l ~=0)    
    A = -F_struc(K.f+1:K.f+K.l,2:end);
    b = F_struc(K.f+1:K.f+K.l,1);        
    n = size(F_struc,2)-1;    
    cand_rows = find(sum(A~=0,2)==1);
    for i = 1:length(cand_rows)
        j = find(A(cand_rows(i),:));
        if A(cand_rows(i),j)>0
            ub(j)=min(ub(j),sup(b(cand_rows(i))/A(cand_rows(i),j)));
        else
            lb(j)=max(lb(j),inf_(b(cand_rows(i))/A(cand_rows(i),j)));
        end
    end
end