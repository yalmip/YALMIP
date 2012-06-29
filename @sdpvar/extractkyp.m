function [A,B,P,M,negated] = extractkyp(sys);
%EXTRACTKYP Returns (A,B,P,M) from KYP object

% Author Johan Löfberg
% $Id: extractkyp.m,v 1.2 2004-07-01 11:17:10 johanl Exp $

A = sys.extra.A;
B = sys.extra.B;
P = sys.extra.P;
M = sys.extra.M;
negated = sys.extra.negated;