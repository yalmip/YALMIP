function varargout = implies(varargin)
%IMPLIES Logical implication
%
% IMPLIES(X,Y) creates a mixed integer representation of
% the constraint X --> Y, i.e. Y is true if X is true.
%
% Syntax
%   F = implies(X,Y,tol)
%
% Input
%   X : binary SDPVAR variable or a constraint
%   Y : binary SDPVAR variable or a constraint
%  tol: Optional threshhold for defining zero (see NOTE)
%
% Output
%   F : SET object
%
% Examples
%
%  binvar X Y; F = implies(X,Y);
%  binvar X;sdpvar Y; F = implies(X,Y>=5);
%  binvar X;Y=sdpvar(3,1); F = implies(X,[sum(Y);Y(2)]>=[5;0]);
%
% Note
%  Using implies with X non-binary is highly sensitive numerically.
%  The problem comes from the definition of 0 in a floating-point
%  environment, and precision in the solver. To account for this,
%  the user can supply a third argument to define a dead-zone around
%  zero, i.e Implies(X<=0,Y) will be replaced with IMPLIES(X<=-tol,Y)
%  Note, you typically need to tweak this number for your
%  application/solver. By default, YALMIP uses tol = 0, which means you
%  easily can get garbage... A positive number means YALMIP is cautious in
%  terms of activating the condition, while a negative number means YALMIP
%  will be aggressive  in activating the condition. 
%
%   See also @SDPVAR/AND, @SDPVAR/OR, IFF

% Author Johan Löfberg
% $Id: implies.m,v 1.6 2009-05-15 10:32:42 joloef Exp $


% There are some cases to take care of...
%
% X --> Y     binary/binary :                     Implemented
% X --> Y     binary/(LP,equality,sdp)            Implemented
% X --> Y     (LP,equality,sdp)/binary            Not implemented
% X --> Y     (LP,equality,sdp)/(LP,equality,sdp) Not implemented

X = varargin{1};
Y = varargin{2};

% % Call recursicely on X -> (A,B,...)
% if isa(varargin{1},'sdpvar') & (isa(varargin{2},'lmi') | isa(varargin{2},'constraint'))
%     if length(varargin{1})==1 & length(varargin{2})>1
%         F = set([]);
%         for i = 1:length(varargin{2})
%             if nargin == 3
%                 F = F + implies(varargin{1},varargin{2}(i),varargin{3});
%             else
%                 F = F + implies(varargin{1},varargin{2}(i));
%             end
%         end
%         varargout{1} = F;
%         return
%     end
% end

if isempty(X)
    varargout{1} = [];
end

switch class(X)

    case {'sdpvar','constraint','lmi'}      
        varargout{1} = setupMeta(lmi([]), mfilename,varargin{:});
        
    case 'char'        
        varargout{1} = implies_internal(varargin{3:end});
        
    case 'logical'
        if length(X)==1
            if X
                varargout{1} = Y;
            else
                varargout{1} = [];
            end
        else
            if length(X) == length(Y)
                i = find(X);
                if isempty(i)
                    varargout{1} = [];
                else
                    varargout{1} = Y(i);
                end
            else
                error('Size mismatch in input arguments');
            end
        end
end

