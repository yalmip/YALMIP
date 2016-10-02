function [merged,map] = mergeDistributions(randomVariables);

merged = {};
map = [];
used = zeros(1,length(randomVariables));
mergedindex = 1;
for i = 1:length(randomVariables)
    if ~used(i)
        this = randomVariables{i};
        used(i) = 1;
        for j = i+1:length(randomVariables)
            if ~used(j)
                that = randomVariables{j};
                if strcmp(func2str(this.distribution.name),'random') && strcmp(func2str(that.distribution.name),'random')
                    if mergable(this.distribution.parameters{1},that.distribution.parameters{1})
                        % Same distribution
                        this = merge2distributions(this,that);
                        map(j) = mergedindex;
                        used(j) = 1;
                    end
                end
            end
        end
        merged{mergedindex} = this;
        map(i) = mergedindex;
        mergedindex = mergedindex + 1;
    end
end
