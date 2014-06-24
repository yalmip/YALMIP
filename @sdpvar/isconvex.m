function YESNONA = isconvex(p)
%ISCONVEX Tries to determine if a scalar expression is convex
%
% T = isconvex(p)
%
% p: scalar SDPVAR object
% T: The result (1: convex, 0: concave, NaN: cannot be determined)
%
% NOTE : Under development. Do not trust if you have monomials in your model
%
% Example
% sdpvar x y
% isconvex(x+y) will return 1
% isconvex(x+y^2) will return 1
% isconvex(exp(x+y)) will return 1
% isconvex(-exp(x+y)) will return 0
% isconvex(max(x,exp(x+y))) will return 1
% isconvex(-max(x,exp(x+y))) will return 0
% isconvex(max(x,min(x,-x))) will return NaN

YESNONA = NaN;
[F,failure,cause] = expandmodel([],p,sdpsettings('allownonconvex',0,'allowmilp',0));
if failure == 0
    YESNONA = true;
else
    [F,failure,cause] = expandmodel([],-p,sdpsettings('allownonconvex',0));
    if failure == 0
        % p is nonconvex
        YESNONA = false;
    end
end
