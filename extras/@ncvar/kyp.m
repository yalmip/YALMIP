function sys = kyp(A,B,P,M)
%KYP Create KYP matrix variable
%
%   X = KYP(A,B,P,M)       
%
% KYP is used to generate the matrix
%
% [A'*P+P*A P*B;B'*P zeros(size(B,2))]+M;
%
% Note, Information is stored internally
% to inform YALMIP that this is a object
% defined from a KYP structure. Hence,
% the objects KYP(A,B,P,M) and the SDPVAR
% [A'*P+P*A P*B;B'*P zeros(size(B,2))]+M
% are not equivalent
%
%  See also @sdpvar/LYAP

% Author Johan Löfberg
% $Id: kyp.m,v 1.1 2006-08-10 18:00:21 joloef Exp $

if isempty(B)
    kyp_part = [A'*P+P*A];
else
    kyp_part = [A'*P+P*A P*B;B'*P zeros(size(B,2))];
end
if isempty(kyp_part)
    sys = M;
else
    if isempty(M)
        sys = kyp_part;
    else
        sys = kyp_part+M;
    end
end
sys.typeflag = 9;
sys.extra.A = A;
sys.extra.B = B;
sys.extra.P = P;
sys.extra.M = M;
sys.extra.negated = 0;