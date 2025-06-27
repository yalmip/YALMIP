function islinear = linearin(p,y)
%LINEARIN Checks if expression is linear w.r.t variables
%
% YESNO = LINEARIN(p,x)
%
% p      : SDPVAR object.
% x      : Linear w.r.t this SDPVAR
%
% Examples
% x1 = sdpvar(1,1);x2 = sdpvar(1,1);
% p = [x1^2 + x1*x2];
%
% linearin(p,x1) returns FALSE
%
% linearin(p,x2) returns TRUE

if isnumeric(p)
    islinear = true;
    return
end

[mt,vt] = yalmip('monomtable');
pv = getvariables(p);
yv = getvariables(y);
mt = mt(pv,:);
degrees = mt(:,yv);
islinear = all(nonzeros(degrees) <= 1);