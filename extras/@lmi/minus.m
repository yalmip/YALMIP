function sys = minus(X,Y)
%MINUS (overloaded)

if isempty(X)
    sys = ([]);
    return
end

if isempty(Y)
    sys = X;
    return
end

if isa(Y,'double') || isa(X,'double')
    error('You cannot substract a point from a set of constraint. Constraint substraction is only used removal of constraints. To translate a set, use REPLACE')
end

idX = getlmiid(X);
idY = getlmiid(Y);
YinX = find(~ismember(idX,idY));
% Get the correct subsref...
sys = subsref(X,struct('type','()','subs',{{YinX}}));
