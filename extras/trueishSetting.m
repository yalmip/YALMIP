function yesno = trueishSetting(q)

if isequal(q,1) || isequal(q,true) || isequal(q,'yes') || iequal(q,'on')
    yesno = true;
else
    yesno = false;
end