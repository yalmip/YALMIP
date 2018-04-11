function varargout=size(varargin)
%SIZE (overloaded)

if nargin == 1        
  bsize  = varargin{1}.dim;
  switch (nargout)
  case 0
    varargout{1} = bsize;
  case 1
    varargout{1} = bsize;
  case 2
    varargout{1} = bsize(1);
    varargout{2} = bsize(2);
  otherwise
    error('>2 outputs in size?');
  end
else
    if varargin{2} > length(varargin{1}.dim)
        varargout{1} = 1;
    else
        varargout{1} = varargin{1}.dim(varargin{2});	
	end
end
