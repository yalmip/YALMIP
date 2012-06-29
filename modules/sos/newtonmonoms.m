function v = newtonmonoms(p)
% NEWTONMONOMS Computes all monoms inside half Newton polytope
%
% V = NEWTONMONOMS(P)
%
% Input
%  P : Scalar SDPVAR object
%
% Output
%  V : Vector with SDPVAR objects
%
% Example:
%
% sdpvar x y
% sdisplay(newtonmonoms(1+x^4*y^2+x^2*y^4))
%
% See also NEWTONREDUCE, CONSISTENT, CONGRUENCEBLOCKS

% Author Johan Löfberg
% $Id: newtonmonoms.m,v 1.1 2006-03-30 13:27:20 joloef Exp $

if isa(p,'double')
    v = 1;
else
    x = recover(depends(p));
    v = monolist(x,degree(p)/2);
    v = newtonreduce(v,p);
end

