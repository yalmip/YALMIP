function sys = mergelmi(varargin)
% MERGELMI            Merges a set of LMI objects
%
% MERGELMI is obsolete and should not be used. Use + instead
%
%    See also   LMI, UPDATELMI, DELLMI

% Author Johan Löfberg
% $Id: mergelmi.m,v 1.3 2005-02-04 10:10:27 johanl Exp $

for i = 1:nargin
    if ~(isa(varargin{i},'lmi'))
        error('All arguments must be LMI objects')
    end
end

sys = varargin{1};
lmitop = size(varargin{1}.clauses,2)+1;
for i = 2:nargin
    F =  varargin{i};
    for j = 1:size(F.clauses,2)
        sys.clauses{lmitop} = F.clauses{j};
        sys.LMIid(lmitop) = F.LMIid(j);
        lmitop = lmitop+1;
    end
end





