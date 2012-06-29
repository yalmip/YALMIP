function F = complements(C1,C2)
%COMPLEMENTS Defines complementary constraints
%   
%   F = COMPLEMENTS(C1,C2)
% 
%   
  
% Author Johan Löfberg
% $Id: complements.m,v 1.4 2009-04-29 12:44:40 joloef Exp $

F = complements(set(C1),set(C2));
	