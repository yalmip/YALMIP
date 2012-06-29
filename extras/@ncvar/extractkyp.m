function [A,B,P,M,negated] = extractkyp(sys);
%EXTRACTKYP Returns (A,B,P,M) from KYP object

% Author Johan Löfberg
% $Id: extractkyp.m,v 1.1 2006-08-10 18:00:20 joloef Exp $

A = sys.extra.A;
B = sys.extra.B;
P = sys.extra.P;
M = sys.extra.M;
negated = sys.extra.negated;