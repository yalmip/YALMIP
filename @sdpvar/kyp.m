function sys = kyp(A,B,P,M)
%KYP Create KYP matrix variable (Legacy)
%
%   X = KYP(A,B,P,M)
%
% KYP is used to generate the matrix
%
% [A'*P+P*A P*B;B'*P zeros(size(B,2))]+M;

if nargin<4
    M = 0;
end

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