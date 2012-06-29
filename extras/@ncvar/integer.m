function x = integer(x)
%INTEGER Constrains variables to be integer
%
%   F = INTEGER(x) is used to posteriori constrain 
%   variables to be integer, in contrast to INTVAR 
%   that declares variables as integer a priori.
%
%   NOTE
%    The integrality constraint is imposed on the involved
%    decision variables, not on the actual SDPVAR object.
%
%   INPUT
%    x : SDPVAR object
%
%   OUTPUT
%    F : SET object
%
%   EXAMPLE
%    F = set(integer(x));             % Full verbose notation
%    F = integer(x);                  % Short notation
%    F = integer(x*pi)                % Equivalent to code above
%    solvesdp(G+set(integer(x)),obj)  % Add integer constraint on the fly
%
%   See also BINARY, SET, SDPVAR, INTVAR, BINVAR

% Author Johan Löfberg
% $Id: integer.m,v 1.1 2006-08-10 18:00:20 joloef Exp $

x.typeflag = 7;
x = lmi(x);