function sys = dellmi(F,thelmi)
%dellmi           (OBSOLETE) Deletes a constraint in an LMI object
%   
%    F = DELLMI(F,n)         Delete the n'th LMIs
%    F = DELLMI(F,'tag')  Delete the LMIs having specified tag
%
%    See also   LMI, UPDATELMI


% Author Johan Löfberg
% $Id: dellmi.m,v 1.3 2005-02-04 10:10:26 johanl Exp $

if nargin ~=2
	error('DELLMI needs two argument')
end

if ~(isa(F,'lmi') & (isa(thelmi,'double') | isa(thelmi,'char')))
	error('First argument should be an lmi object and second argument integer or string')
end

% If indexed using handle, convert to index
if isa(thelmi,'char')
	thelmitemp = [];
	for i = 1:size(F.clauses,2)
		if strcmp(F.clauses{i}.handle,thelmi)
			thelmitemp=[thelmitemp i];
		end
	end
	if isempty(thelmitemp)
		em = ['LMI ''' thelmi ''' not available.'];
		error(em)
	else
		thelmi = thelmitemp;
	end
end

% Checks so that it exist
if any((thelmi<1)) | any((thelmi>size(F.clauses,2)))
	em = ['LMI #' num2str(thelmi) ' not available.'];
	error(em)
end
sys = F;
del_lmi = thelmi(logical(thelmi<=size(F.clauses,2)));
if ~isempty(del_lmi)
	keep_lmi = setdiff(1:size(F.clauses,2),del_lmi);
	sys.clauses = {F.clauses{keep_lmi}};
    sys.LMIid = sys.LMIid(keep_lmi);
end
