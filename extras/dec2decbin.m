function s=dec2bin(d,n)
%DEC2BIN Internal function generate binary matrices

% Author Johan Löfberg
% $Id: dec2decbin.m,v 1.2 2004-07-02 08:17:30 johanl Exp $

[f,e]=log2(max(d)); % How many digits do we need to represent the numbers?
s=rem(floor(d(:)*pow2(1-max(n,e):0)),2);
