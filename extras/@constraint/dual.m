function sys = dual(X)
%DUAL Extract dual variable
%   
%   Z = DUAL(F)     Returns the dual variable for the constraint F
% 
%   See also SOLVESDP, DUALIZE
  
% Author Johan Löfberg
% $Id: dual.m,v 1.2 2010-02-08 13:06:11 joloef Exp $

sys = dual(set(X));
