function x = binary(x)
%BINARY Constrains variables to be binary (0/1)
%
%   F = BINARY(x) is used to a posteriori constrain 
%   variables to be binary, in contrast to BINVAR 
%   that declares variables as binary a priori.
%
%   NOTE
%    The binary constraint is imposed on the involved
%    decision variables, not on the actual SDPVAR object.
%
%   INPUT
%    x : SDPVAR object
%
%   OUTPUT
%    F : SET object
%
%   EXAMPLE
%    F = binary(x);                  
%    F = binary(x*pi)   % Equivalent to code above
%
%   See also INTEGER, SDPVAR, INTVAR, BINVAR

x.typeflag = 8;
x = lmi(x);