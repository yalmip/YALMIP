function self = sample(self,N)
%SAMPLE Draw a sample in an optimizer object and instantiates parameter
%
%   Q = SAMPLE(P,N) generates concatenated instantiated optimizer objects
%   where samples have been drawn for all parameters with associated
%   samplers.
%
%   INPUT
%    P : OPTIMIZER object
%    N : Number of samples to draw
%
%   OUTPUT
%    Q : OPTIMIZER object
%
%   EXAMPLE
%    sdpvar x(2,1) w(2,1)
%    F = [-10 <= x + w <= 10, uncertain(w,'unif',2,3)]
%    P = optimizer(F,sum(x),[],w,x)
%    Q = sample(P,10);
%    plot(Q);
%    xoptimal = Q([]);
%
%   See also OPTIMIZER, UNCERTAIN

if nargin < 2
    N = 1;
end
allSamples = {};
for k = 1:N
    cells = cell(1,length(self.diminOrig));
    for i = 1:length(self.diminOrig)
        if isfield(self.input,'stochastics')
            if ~isempty(self.input.stochastics{i});
                temp = {self.input.stochastics{i}.name,self.input.stochastics{i}.parameters{:},self.diminOrig{i}};
                sampledData = feval(temp{:});
                cells{i} = sampledData;
            else
                cells{i} = [];
            end
        end
    end
    if ~isempty(cells)
        y.type = '{}';
        y.subs = {cells{:},'nosolve'};      
        allSamples{end + 1} = subsref(self,y);
    else
        return
    end
end
self = horzcat(allSamples{:});


