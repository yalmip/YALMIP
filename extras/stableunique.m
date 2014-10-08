function [C,firstPos,index,skipped] = stableunique(varargin)

% Implements unique(a,'rows','stable') for a vector
% Needed to support older versions of MATLAB
% Brute-force implementation, no brain involved.

a = varargin{1};
n = length(a);
% Simple cases
if n == 1
    C = a;
    firstPos=1;index=1;
    skipped = [];
elseif all(diff(sort(a))>0)
    C = a;
    firstPos = 1:n;
    index = 1:n;
    skipped = [];
end
C = a(1);
firstPos = 1;
index = zeros(length(a),1);
index(1) = 1;
skipped = [];
for i = 2:length(a)
    pos = min(find(a(i)==C));
    if isempty(pos)
        C = [C;a(i)];
        firstPos = [firstPos;i];
        pos = length(C);
    else
        skipped = [skipped i];
    end
    index(i)=pos;
end
