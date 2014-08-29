function s=dec2decbin(d,n)
%DEC2BIN Internal function generate binary matrices

[f,e]=log2(max(d)); % How many digits do we need to represent the numbers?
s=rem(floor(d(:)*pow2(1-max(n,e):0)),2);
