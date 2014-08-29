function [KKTConstraints, details] = kkt(P)
%DUAL Create KKT system for optimizer P
%
% [KKTConstraints, details] = kkt(P)

[KKTConstraints, details] = kkt(P.F,P.h,P.output.expression);

