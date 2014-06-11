function varargout=size(varargin)
%size              Returns the number of inequalities and equalities 
%   
%    n = SIZE(F)     Returns the number of inequalities
%    [n,m] = SIZE(F) Returns the number of inequalities and equalities
%
%    See also   LMI, ADDLMI

F = varargin{1};
nequ = 0;
nlmi = 0;
for i = 1:size(F.clauses,2)
    if F.clauses{i}.type==3
        nequ = nequ + 1;
    else
        nlmi = nlmi + 1;
    end
end

if nargin == 1  
  switch (nargout)
  case {0,1}
    varargout{1} = nlmi;
  case 2
    varargout{1} = nlmi;
    varargout{2} = nequ;
  otherwise
    error('>2 outputs in size?');
  end
else
	switch(varargin{2})
	case 1
		varargout{1} = nlmi;
	case 2
		varargout{1} = nequ;
	otherwise
		error('Second argument should be 1 (# LMIs) or 2 (# equalities)');
	end
end