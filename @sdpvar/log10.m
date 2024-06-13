function varargout = log10(varargin)

varargout{1} = log(varargin{1})*(1/log(10));