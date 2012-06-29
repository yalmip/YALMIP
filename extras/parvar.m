function sys = parvar(varargin)
%PARVAR Create symbolic parametric variable 
%   
%   PAPRVAR works exactly as SDPVAR, with the only difference that
%   the elements in the variable automatically will be defined
%   as parametric variables (for a multi-parametric program)
%
%   See also PARAMETRIC, INTVAR, BINVAR, SDPVAR, BINARY, INTEGER

% Author Johan Löfberg
% $Id: parvar.m,v 1.2 2004-07-09 14:29:39 johanl Exp $

switch nargin
    case 1
        sys = sdpvar(varargin{1});
    case 2
        sys = sdpvar(varargin{1},varargin{2});
    case 3
        sys = sdpvar(varargin{1},varargin{2},varargin{3});
    case 4
        sys = sdpvar(varargin{1},varargin{2},varargin{3},varargin{4});
    otherwise
        error('Wrong number of input arguments. See help-text for sdpvar')
end

yalmip('setparvariables',[yalmip('parvariables') getvariables(sys)]);