function varargout = not(varargin)
%NOT (overloaded)
%   
%    y = not(x)       
%    y = ~x    
%
% Short for y = 1-x.
%
% It is assumed that x is a binary variable (either explicitely declared
% using BINVAR, or constrained using BINARY.)
%
%    See also SDPVAR/OR, SDPVAR/AND, BINVAR, BINARY

% Author Johan Löfberg 
% $Id: not.m,v 1.6 2005-02-10 12:26:38 johanl Exp $   

varargout{1}=1-varargin{1};
