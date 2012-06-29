function Z = minus(P)
%display           Overloaded

% Author Johan Löfberg 
% $Id: uminus.m,v 1.1 2004-06-17 08:40:09 johanl Exp $  

Z = P;
Z.cx =  - Z.cx;
Z.gain = -Z.gain;
