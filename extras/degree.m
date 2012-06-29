function deg = degree(x,y,z)
%DEGREE Polynomial degree
%
% DEG = DEGREE(p,x,e)
%
% p : SDPVAR object
% x : Degree w.r.t linear SDPVAR objects, can be [].

% e : If e=1, returns degree of each element in p
%
% Examples
% x1 = sdpvar(1,1);x2 = sdpvar(1,1);
% p = [x1;x1*x2+x2^2];
%
% degree(p) returns 2
%
% degree(p,x1) returns 1
%
% degree(p,[x1 x2]) returns [1 2]
%
% degree(p,[x1 x2],1) returns [1 0;1 2]
%
% degree(p,[],1) returns [1;3]

% Author Johan Löfberg 
% $Id: degree.m,v 1.2 2004-07-02 08:17:30 johanl Exp $  

if isa(x,'double')
    deg = 0;
else
    error('DEGREE only defined for SDPVAR and DOUBLE.');
end
    