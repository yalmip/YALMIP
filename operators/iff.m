function varargout = iff(varargin)
%IFF Logical equivalence
%
% IFF(X,Y) creates a mixed integer representation of
% the constraint X <--> Y, i.e. Y is true iff X is true.
%
% Syntax
%   F = iff(X,Y)
%
% Input
%   X : binary SDPVAR variable or a set of linear (in)equalities
%   Y : binary SDPVAR variable or a set of linear (in)equalities
%
% Output
%   F : Constraint object
%
% Examples
%
%  binvar X,Y; F = set(iff(X,Y));
%  sdpvar X;binvar Y; F = set(iff(X>5,Y));
%  sdpvar X;binvar Y; F = set(iff(Y,X==5));
%
% Overloading
%
% The iff overloads == for logic constraints.
%
%  sdpvar X;binvar Y; F = set((X>=5) == Y);
%  sdpvar X;binvar Y; F = set(Y == (X==5));

%
% Note
%  The function IFF is not complete, but will be
%  improved upon in future releases.
%
%   See also @SDPVAR/AND, @SDPVAR/OR, IMPLIES

% Author Johan Löfberg
% $Id: iff.m,v 1.4 2007-08-02 19:17:36 joloef Exp $

X = varargin{1};
Y = varargin{2};

switch class(varargin{1})
    case {'lmi','constraint','sdpvar'}
        varargout{1} = setupMeta(lmi([]), mfilename,varargin{:});


    case 'char'
        varargout{1} = iff_internal(varargin{3:end});
end