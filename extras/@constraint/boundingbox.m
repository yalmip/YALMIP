function varargout = boundingbox(F,ops,x)
%BOUNDINGBOX Computes bounding box of a constraint object
%
% If three outputs are requrested, the numerical bounds are returned
% [C,L,U] = boundingbox(F)
%
% If only one output is requested, only the symbolic model is requested
% C = boundingbox(F)

% $Id: boundingbox.m,v 1.2 2008-02-14 14:53:36 joloef Exp $

if nargin < 2
    ops = sdpsettings('verbose',0);
end

if nargin < 3
    [B,L,U] = boundingbox(set(F),ops);
else
    [B,L,U] = boundingbox(set(F),ops,x);
end
switch nargout
    case 0
        B
    case 1
        varargout{1} = B;
    case 2
        varargout{1} = B;
        varargout{2} = L;
    case 3
        varargout{1} = B;
        varargout{2} = L;
        varargout{3} = U;
    otherwise
end
