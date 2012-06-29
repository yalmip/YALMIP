function y = rdivide(X,Y)
%RDIVIDE (overloaded)

% Author Johan Löfberg 
% $Id: rdivide.m,v 1.7 2006-01-26 13:44:13 joloef Exp $   

% Check dimensions
[nx,mx] = size(X);
[ny,my] = size(Y);
if ~((prod(size(X))==1) | (prod(size(Y))==1))
    if ~((nx==ny & (mx ==my)))
        error('Matrix dimensions must agree.')
    end
end

% Quick exit for simple case X/scalar
if isa(Y,'double') & prod(size(Y))==1
    y = X;
    y.basis = y.basis/Y;
    % Reset info about conic terms
    y.conicinfo = [0 0];
    return
end

if isa(X,'sdpvar') & isa(Y,'double')
    y = X.*(1./Y);
    return
end

% FIX : SLOOOOW BUT SOMEWHAT ROBUST
[nx,mx] = size(X);
[ny,my] = size(Y);
if prod(size(X)) == 1 & prod(size(Y))~=1
    X = repmat(X,ny,my);
end;
if prod(size(Y)) == 1 & prod(size(X))~=1
    Y = repmat(Y,nx,mx);
end;
[nx,mx] = size(X);
y = [];
for i = 1:nx   
    if mx==1
        dummy = struct('type','()','subs',{{i,1}});
        y=[y;subsref(X,dummy)*subsref(Y,dummy)^-1];
    else
        ytemp = [];
        for j = 1:mx
            dummy = struct('type','()','subs',{{i,j}});
            ytemp = [ytemp subsref(X,dummy)*subsref(Y,dummy)^-1];
        end
        y = [y;ytemp];
    end
end
% Reset info about conic terms
if isa(y,'sdpvar')
    y.conicinfo = [0 0];
end