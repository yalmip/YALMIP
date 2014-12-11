function H = deriveBasis(A_equ)
[L,U,P] = lu(A_equ);
[L,U,P] = lu(A_equ');
r = colspaces(L');
AA = L';
H1 = AA(:,r);
H2 = AA(:,setdiff(1:size(AA,2),r));
H = P'*[-H1\H2;speye(size(H2,2))];
function  [indx]=colspaces(A)
indx = [];
for i = 1:size(A,2)
    s = max(find(A(:,i)));
    indx = [indx s];
end
indx = unique(indx);  