function [KKTConstraints, details] = kkt(P)
%DUAL Create KKT system for optimizer P
%
% [KKTConstraints, details] = kkt(P)

% Author Johan Löfberg
% $Id: kkt.m,v 1.2 2010-02-08 13:06:11 joloef Exp $

[KKTConstraints, details] = kkt(P.F,P.h,P.output.expression);

