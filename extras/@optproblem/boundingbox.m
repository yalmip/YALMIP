function [B,L,U] = boundingbox(P,Options,x)
%BOUNDINGBOX Computes bounding box of a constraint object
%
% If only one output is requested, only the symbolic model is returned
% B = boundingbox(F)
%
% If three outputs are requrested, the symbolic model is appended
% [B,L,U] = boundingbox(F)
%
% A second argument can be used to specificy solver settings
% B = boundingbox(F,sdpsettings('solver','cplex'))
%
% A third argument can be used to obtain a bounding box in a particular 
% set of variables, i.e. the bounding box of a projection
% B = boundingbox(F,[],x)

if nargin < 2
    Options = sdpsettings('verbose',0);
end

if nargin == 3
    [B,L,U] = boundingbox(set(P.Constraints),Options,x)
else
    [B,L,U] = boundingbox(set(P.Constraints),Options)
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