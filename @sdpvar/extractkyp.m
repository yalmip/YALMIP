function [A,B,P,M,negated] = extractkyp(sys);
%EXTRACTKYP Returns (A,B,P,M) from KYP object

A = sys.extra.A;
B = sys.extra.B;
P = sys.extra.P;
M = sys.extra.M;
negated = sys.extra.negated;