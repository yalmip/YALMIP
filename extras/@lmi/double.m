function sys = double(X)
%double           Overloaded.

% Author Johan Löfberg
% $Id: double.m,v 1.3 2005-02-12 13:24:38 johanl Exp $

nlmi = size(X.clauses,2);

if (nlmi == 0) 
 sys = NaN;
end

if nlmi>1 
 error('Double not applicable on list of constraints')
end

sys = double(X.clauses{1}.data);
