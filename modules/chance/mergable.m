function yesno = mergable(a,b)
if isequal(a,b)
    yesno = 1;
elseif any(strcmp(a,{'normal','normalf','normalm'})) && any(strcmp(b,{'normal','normalf','normalm'}))
    yesno = 1;
else
    yesno = 0;
end