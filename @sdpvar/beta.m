function varargout=beta(varargin)

if isa(varargin{1},'sdpvar') && isa(varargin{2},'sdpvar')
    error('W and Z can not both be an SDPVAR objects')
end
if isa(varargin{1},'sdpvar')
    varargout{1} = InstantiateElementWise('beta_z',varargin{:});
else
    varargout{1} = InstantiateElementWise('beta_w',varargin{2},varargin{1});
end
