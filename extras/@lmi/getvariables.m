function used = getvariables(F)

F = flatten(F);
if length(F.clauses) == 0
    used = [];
    return
end

if isa(F.clauses{1},'cell')
    F = flatten(F);
end

used = recursivegetvariables(F,1,length(F.clauses));

function used = recursivegetvariables(F,startindex,endindex)

if endindex-startindex>50
    newstart = startindex;
    mid = ceil((startindex + endindex)/2);
    newend = endindex;
    used1 = recursivegetvariables(F,newstart,mid);
    used2 = recursivegetvariables(F,mid+1,newend);
    used = uniquestripped([used1 used2]);
else
    used = [];
    if startindex <= length(F.clauses)
        used = getvariables(F.clauses{startindex}.data);
        for i = startindex+1:endindex
            Fivars = getvariables(F.clauses{i}.data);
            if ~isequal(used,Fivars(:)')
                used = [used Fivars(:)'];
            end
        end
        used = uniquestripped(used);
    end
end