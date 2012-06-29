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
%    F = set(binary(x));             % Full verbose notation
%    F = binary(x);                  % Short notation
%    F = binary(x*pi)                % Equivalent to code above
%    solvesdp(G+set(binary(x)),obj)  % Add binary constraint on the fly
%
%   See also INTEGER, SET, SDPVAR, INTVAR, BINVAR

% Author Johan Löfberg
% $Id: binary.m,v 1.1 2006-08-10 18:00:19 joloef Exp $

x.typeflag = 8;
x = lmi(x);