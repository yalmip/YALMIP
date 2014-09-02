function C = tplus(varargin)
% TPLUS Tropical addition (max-plus algebra)
%
% C = TPLUS(A,B) Computes tropical addition C(i,j) = max(A(i,:) + B(:,j)');
%
% Note that TPLUS can just as well can be replaced with C = MAX(A,B)
%
% See also ttimes

if nargin == 2
C = max(varargin{1},varargin{2});
else
    C = max(varargin{1},varargin{2});
    for i = 3:nargin
        C = max(C,varargin{i});
    end
end
