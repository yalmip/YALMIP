function varargout = plothyperplanes(varargin)
%PLOTLATTICE  Draws the hyperplanes of a polytope
%
% p = plothyperplanes(P,color,D)
%
% Note that only polytopes in R^2 are supported.
%
% C    :  Constraint object
% color:  color [double] ([r g b] format) or char from 'rymcgbk'
% D    :  Region over which plot is generated. Based on bounding box if not supplied

P = varargin{1};

if length(depends(P)) > 2
    error('PLOTHYPERPLANES is only applicable in 2D');
end

if nargin > 1
    color = varargin{2};
    if isempty(color)
        color = 'blue';
    end
else
    color = 'blue';
end
if nargin > 2
    domain = varargin{3};
    X = recover(depends(domain));
    [~,L,U] = boundingbox(domain);
    oldL = L;
    oldU = U;
else
    domain = P;
    X = recover(depends(P));
    [~,L,U] = boundingbox(P);
    oldL = L;
    oldU = U;
    % Explode the box
    C = (U+L)/2;
    W = (U-L)/2;
    U = C+1.5*W;
    L = C-1.5*W;       
end

bA = getbase(sdpvar(P));
ls = [];
for i = 1:size(bA,1);
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
    ls = [ls l];
end    
% L = oldL;
% U = oldU;
% C = (L+U)/2;W = (U-L)/2;
% L = C-2*W;
% U = C+2*W;
axis([L(1) U(1) L(2) U(2)])
varargout{1} = ls;
    
    
    