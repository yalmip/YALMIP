function F = replace(F,X,W,expand)
%REPLACE        Substitutes variables in a SET object
%
%Z = REPLACE(F,X,W)  Replaces any occurence of the SDPVAR object Y 
%                    in the SDPVAR object X with the double W
%
% Example
%  x = sdpvar(1,1);
%  t = sdpvar(1,1);
%  F = set([1+t;1+x+t] > 0);
%  F = replace(F,x,2) generates F=set([1+t;3+t] > 0)

if nargin<4
    expand = 1;
end

if prod(size(W)) == 1
    W = repmat(W,size(X));
end

keep = [];

for i = 1:length(F.clauses)
    F.clauses{i}.data = replace(F.clauses{i}.data,X,W,expand);
    if isempty(F.clauses{i}.data) | isa(F.clauses{i}.data,'double')
        keep(i)=0;
    else
        keep(i)=1;
    end
end

F.clauses = {F.clauses{find(keep)}};

% Get new identifiers
F.LMIid = [];
for i = 1:length(F.clauses)
    F.LMIid = [F.LMIid yalmip('lmiid')];
end