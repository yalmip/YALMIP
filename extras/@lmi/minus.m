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

idX = getlmiid(X);
idY = getlmiid(Y);
YinX = find(~ismember(idX,idY));
% Get the correct subsref...
sys = subsref(X,struct('type','()','subs',{{YinX}}));
