function F = horzcat(varargin)

if all(cellfun('isclass',varargin,'lmi'))
    F = fastcat(varargin{:});
    return
end

F = [];
for i=1:1:nargin
    if isa(varargin{i},'double') & ~isempty(varargin{i})
        warning('One of the constraints evaluates to a DOUBLE variable');
    elseif isa(varargin{i},'logical')
        if all(varargin{i}==1)
            %  warning('One of the constraints evaluates to a LOGICAL variable');
        else
            error('One of the constraints evaluates to a FALSE LOGICAL variable');
        end
    elseif isa(varargin{i},'optproblem')
        F = [varargin{i},F,varargin{i+1:end}];
        return;
    else
        H = set(varargin{i});
        F = F + H;
    end
end

function X = fastcat(varargin)

X = varargin{1};
nTOT = length(X.clauses);
for i = 2:nargin
    X.clauses = cat(2,X.clauses,varargin{i}.clauses);
    nTOT = nTOT + length(varargin{i}.clauses);
    X.LMIid = [X.LMIid varargin{i}.LMIid];
end

% VERY FAST UNIQUE BECAUSE THIS IS CALLED A LOT OF TIMES....
i = sort(X.LMIid);
i = i(diff([i NaN])~=0);
if length(i)<nTOT
    [i,j] = unique(X.LMIid);
    X = subsref(X,struct('type','()','subs',{{j}}));
end
