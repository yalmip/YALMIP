function self = sample(self,N)
if nargin < 2
    N = 1;
end
allSamples = {};
for k = 1:N
    cells = cell(1,length(self.diminOrig));
    for i = 1:length(self.diminOrig)
        if ~isempty(self.input.stochastics{i});
            temp = {@random,self.input.stochastics{i}.name,self.input.stochastics{1}.parameters{:},self.diminOrig{i}};
            sampledData = feval(temp{:});
            cells{i} = sampledData;
        else
            cells{i} = [];
        end
    end
    y.type = '{}';
    y.subs{1} = cells;
    y.subs{2} = 'nosolve';
    allSamples{end + 1} = subsref(self,y);
end
self = horzcat(allSamples{:});


