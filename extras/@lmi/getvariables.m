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
        
        if F.clauses{startindex}.type == 56
            used = [];
            for j = 1:length(F.clauses{startindex}.data)
                used = [used getvariables(F.clauses{startindex}.data{j})];
            end
        else
            used = getvariables(F.clauses{startindex}.data);
        end
        
        for i = startindex+1:endindex
            
            %Fivars = getvariables(F.clauses{i}.data);
            if F.clauses{i}.type == 56
                % Meta constraint such as implies. This object is just holding
                % the data involved
                Fivars = [];
                for j = 1:length(F.clauses{i}.data)
                    Fivars = [Fivars getvariables(F.clauses{i}.data{j})];
                end
            else
                Fivars = getvariables(F.clauses{i}.data);
            end
                        
            if ~isequal(used,Fivars(:)')
                used = [used Fivars(:)'];
            end
        end
        used = uniquestripped(used);
    end
end