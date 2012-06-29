function Y = uminus(Y)
%UMINUS (overloaded)

% Author Johan Löfberg
% $Id: uminus.m,v 1.7 2009-08-26 01:46:41 joloef Exp $

if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

switch Y.typeflag
    case {0,1,2,3,4}
        Y.basis = -Y.basis;
    case {9,40} % Simple KYP, to be obsoleted       
        Y.basis = -Y.basis
        Y.extra.M = -Y.extra.M;
        Y.extra.negated = ~Y.extra.negated;
%     case 40 % Generalized KYP
%         Y.basis = -Y.basis;
%         Y.extra.M = -Y.extra.M;
%         for i = 1:length(Y.extra.K)
%             Y.extra.negated(i) = ~Y.extra.negated(i);
%         end
    case 5
        error('Cone object cannot be negated');
    otherwise
end
% Reset info about conic terms
Y.conicinfo = [0 0];
Y.extra.opname='';
Y = negatefactors(Y);


