function sys = blocks(varargin)
% BLOCKS Internal function (Not used)
%
% BLOCKS(X1,X2,X3,X4)  Tries to create a blocks-symmetric matrix
%
% Example: blocks(X1,X2,X3) creates [X1 X2;X' X3]
%          blocks(X1,X2,X3,X4,X5,X6) creates [X1 X2 X3;X2' X4 X5;X3' X5' X6]

if rem(nargin,3)~=0
    error('The number of blocks in a block-symmetric matrix has to be a multiple of 3')
end

m = nargin;
n = max(roots([1 1 -2*m]));

k = 1;
Z = zeros(n);
for i = 1:n
    for j = i:n
        Z(i,j)=k;k = k+1;
    end
end
Z = (Z+Z')-diag(diag(Z));

sys = [];
for i = 1:n
    sysrow = [];
    for j = 1:1:n  
        if j<i
            sysrow = [sysrow varargin{Z(i,j)}];
        else
            sysrow = [sysrow varargin{Z(i,j)}'];     
        end
    end
    sys = [sys;sysrow];
end