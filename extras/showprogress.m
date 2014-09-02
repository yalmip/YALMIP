function showprogress(thetext,doit)
%SHOWPROGRESS Internal function for printing messages

if doit>0
	fprintf('+ %s\n',thetext);
end
