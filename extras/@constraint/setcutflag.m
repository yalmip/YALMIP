function X = setcutflag(X,flag)
%setcutlag Internal : defines a SET object as a CUT

X = flatten(X);
if nargin == 1
    flag = 1;
end
for i = 1:length(X.clauses)
    X.clauses{i}.cut = flag;
end