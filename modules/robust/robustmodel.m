function varargout = robustmodel(varargin)
%ROBUSTMODEL  Derives robust counterpart.
%
% [Frobust,objrobust,failure] = ROBUSTMODEL(F,h,options) is used to derive
% the robust counterpart of an uncertain YALMIP model.
%
%   min        h(x,w)
%   subject to
%           F(x,w) >(=) 0  for all w in W
%
% The constraints and objective have to satisfy a number of conditions for
% the robustification to be possible. Please refer to the YALMIP Wiki for
% the current assumptions.
%
% Some options for the robustification strategies can be altered via the
% solver tag 'robust' in sdpsettings
%
%  'robust.lplp'  : Controls how linear constraints with affine
%                   parameterization in an uncertainty with polytopic
%                   description is handled. Can be either 'duality' or
%                   'enumeration' 
%
%  'robust.auxred': Controls how uncertainty dependent auxiliary variables
%                   are handled
%                   Can be either 'projection' or 'enumeration' (exact),
%                   or 'none' or 'affine' (conservative)
%
%  'robust.reducedual' Controls if the system equality constraints derived
%                   when using the duality filter should be eliminated,
%                   thus reducing the number of variables, possibly
%                   destroying sparsity .
%
%  'robust.polya'  : Controls the relaxation order of polynomials. If set to
%                   NAN, the polynomials will be eliminated by forcing the
%                   coefficients to zero
%
% See also UNCERTAIN

[varargout{1:nargout}] = robustify(varargin{:});