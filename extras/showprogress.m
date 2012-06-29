function showprogress(thetext,doit)
%SHOWPROGRESS Internal function for printing messages

% Author Johan Löfberg
% $Id: showprogress.m,v 1.2 2004-07-02 08:17:32 johanl Exp $

if doit>0
	fprintf('+ %s\n',thetext);
end
