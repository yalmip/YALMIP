function varargout=size(varargin)
%SIZE (overloaded)

% Author Johan Löfberg 
% $Id: size.m,v 1.1 2006-08-10 18:00:22 joloef Exp $   

if nargin == 1    
  bsize  = varargin{1}.dim;%[varargin{1}.n varargin{1}.m];
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
	switch varargin{2}
	case 1
		varargout{1} = varargin{1}.dim(1);
	case 2
		varargout{1} = varargin{1}.dim(2);
	otherwise
		error('Report bug in size')
	end
end
