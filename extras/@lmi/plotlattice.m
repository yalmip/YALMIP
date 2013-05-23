function varargout = plotlattice(varargin)
%PLOTLATTICE  Plots an integer lattice
%
% p = plotlattice(C,which,c,size,options)
%
% Note that only convex sets C in R^2 are supported.
%
% C    :  Constraint object
% which:  'inner' or 'outer'
% color:  color [double] ([r g b] format) or char from 'rymcgbk'
% size :  Size of marker
% options: options structure from sdpsettings
% Example
% sdpvar x1 x2
% plot(x1^2+x2^2 <= 1.5,'green');hold on
% plotlattice(x1^2+x2^2 <= 1.5,'outer','yellow');
% plotlattice(x1^2+x2^2 <= 1.5,'inner','black');

F = varargin{1};
if nargin > 1
    which = varargin{2};
else
    which = 'outer';
end
if nargin > 2
    color = varargin{3};
else
    color = 'yellow';
end
if nargin > 3
    size = varargin{4};
else
    size = 5;
end
if nargin > 4    
    ops = varargin{5};
    if ~isempty(ops)
        ops = sdpsettings(ops,'relax',2,'verbose',0);
    else
        ops = sdpsettings('relax',2,'verbose',0);    
    end
else
    ops = sdpsettings('relax',2,'verbose',0);
end

X = recover(depends(F));
[~,L,U] = boundingbox(F,ops);

x = floor(L(1)):ceil(U(1));
y = floor(L(2)):ceil(U(2));

for i = x
    for j = y
        switch which
            case 'outer'
                l = plot(i,j,'or','MarkerSize',size);
                set(l,'MarkerFaceColor',color);
            case 'inner'
                assign(X,[i;j]);
                p = checkset(F);
                if min(p) >= 0
                    l = plot(i,j,'or','MarkerSize',size,'MarkerFaceColor','yellow')
                    set(l,'MarkerFaceColor',color);
                end
            otherwise
                error
        end
    end
end