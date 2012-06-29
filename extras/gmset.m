function F = gmset(t,x1,x2)
%GMSET  Internal function used for MAXDET formulation

% Author Johan Löfberg
% $Id: gmset.m,v 1.2 2004-07-02 08:17:31 johanl Exp $


F =  set(cone([t;(x1-x2)/2],(x1+x2)/2));