function p = depends(x)
%DEPENDS Returns indicies to variables used in an SDPVAR object
%
% i = depends(x)
%
% Input
%    x : SDPVAR object
% Output
%    i : DOUBLE

% Author Johan Löfberg 
% $Id: depends.m,v 1.2 2004-07-02 08:17:30 johanl Exp $  

p=[];
