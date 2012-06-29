function C = tplus(varargin)
% TPLUS Tropical addition (max-plus algebra)
%
% C = TPLUS(A,B) Computes tropical addition C(i,j) = max(A(i,:) + B(:,j)');
%
% Note that TPLUS can just as well can be replaced with C = MAX(A,B)
%
% See also ttimes

% Author Johan Löfberg 
% $Id: tplus.m,v 1.1 2006-11-24 13:39:43 joloef Exp $  

if nargin == 2
C = max(varargin{1},varargin{2});
else
    C = max(varargin{1},varargin{2});
    for i = 3:nargin
        C = max(C,varargin{i});
    end
end
% [];
% for i = 1:size(A,1)
%     temp = [];
%     for j = 1:size(B,2)
%         temp = [temp max(A(i,j),B(i,j))];
%     end
%     C = [C;temp];
% end
