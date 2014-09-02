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
%    F : Constraint object
%
%   EXAMPLE
%    F = integer(x);                  % Add integrality
%    F = integer(x*pi)                % Equivalent to code above
%
%   See also BINARY, SET, SDPVAR, INTVAR, BINVAR

x.typeflag = 7;
x = lmi(x);