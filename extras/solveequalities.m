function [x_equ,H,A_equ,b_equ,factors] = solveequalities(F_struc,K,unitary)
%SOLVEEQUALITIES Internal function remove equality constraints

% Extract the inequalities
A_equ = F_struc(1:K.f,2:end);
b_equ = -F_struc(1:K.f,1);

if nargin<3
    unitary = 1;
end

factors = [];

if ~unitary

    % This code needs a clean up, and I need to take a course in numerical
    % analysis again.

    % Remove redundant constraints
      [L,U,P] = lu(A_equ);

    % Find a basis for the column space of A_equ
    [L,U,P] = lu(A_equ');
    %  U'*L'*P = A_equ
    %  U'*L'*Px = b_equ
    %  L'*Px = (U')\b_equ
    %  L'*z = (U')\b_equ, z = Px
    %  [FULLRANK other]z = b
    r = colspaces(L');
    AA = L';
    H1 = AA(:,r);
    H2 = AA(:,setdiff(1:size(AA,2),r));
    try
        %x_equ = A_equ\b_equ;
        x_equ = P'*linsolve(full(L'),linsolve(full(U'),full(b_equ),struct('LT',1==1)),struct('UT',1==1));
    catch
        x_equ = A_equ\b_equ;
    end
    % FIX : use L and U stupid!
    H = P'*[-H1\H2;speye(size(H2,2))];
else
    %Improve numerics by removing obviously redundant constraints. This
    %speeds up some cases too
    hash = 1+rand(1+size(A_equ,2),1);
    [i,j] = unique([A_equ b_equ]*hash);
    remove =  setdiff(1:size(A_equ,1),j);
    if ~isempty(remove)
        A_equ(remove,:) = [];
        b_equ(remove) = [];
    end
    % Use unitary basis
    try
        [Q,R,E] = qr(full(A_equ)');
    catch
        [Q,R,E] = qr(A_equ'); % Ouch, that big!
    end
    n = max(find(sum(abs(R),2)>1e-14*size(R,2)));

    Q1 = Q(:,1:n);
    R = R(1:n,:);
    x_equ = Q1*(R'\E'*b_equ);
    H = Q(:,n+1:end); % New basis
end

function  [indx]=colspaces(A)
indx = [];
for i = 1:size(A,2)
    s = max(find(A(:,i)));
    indx = [indx s];
end
indx = unique(indx);