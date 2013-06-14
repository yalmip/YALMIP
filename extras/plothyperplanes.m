function varargout = plothyperplanes(varargin)
%PLOTLATTICE  Draws the hyperplanes of a polytope
%
% p = plothyperplanes(P,color)
%
% Note that only polytopes in R^2 are supported.
%
% C    :  Constraint object
% color:  color [double] ([r g b] format) or char from 'rymcgbk'

F = varargin{1};
if nargin > 1
    color = varargin{2};
else
    color = 'blue';
end
 
P = varargin{1};
X = recover(depends(P));
[~,L,U] = boundingbox(P);

% Explode the box
C = (U+L)/2;
W = (U-L)/2;
U = C+5*W;
L = C-5*W;

bA = getbase(sdpvar(P));
for i = 1:length(bA);
    b = bA(i,1);
    a = -bA(i,2:end);
    if a(1) == 0
        p1 = [L(1) b/a(2)];
        p2 = [U(1) b/a(2)];
    elseif a(2) == 0
        p1 = [b/a(1) L(2)];
        p2 = [b/a(1) U(2)]; 
    else
        p1 = [L(1) (b-a(1)*L(1))/a(2)];
        p2 = [U(1) (b-a(1)*U(1))/a(2)];
    end
    l = line([p1(1) p2(1)],[p1(2) p2(2)]);
    set(l,'color',color);
end    
    
    
    