function F = rational(x,N)
%INTEGER Constrains variables to be rational
%
%   F = RATIONAL(x,N) is used to constrain a variable to be rational in
%   the form x = y/N where y is integer
%
%   For example, rational(x,10) will allow us to search over the numbers
%   0,+-1/10, +-2/10, +- 3/10,...

%   See also BINARY, INTEGER, BINVAR, INTVAR, SEMIVAR

y = intvar(numel(x),1);
F = [reshape(x,[],1) == y/N];