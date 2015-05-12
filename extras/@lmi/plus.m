function X = plus(X,Y)
%PLUS Merges two LMI objects to one LMI

if isa(X,'constraint')
    X = lmi(X);
elseif isa(X,'sdpvar')
    X = lmi(X);
end

if isa(Y,'constraint')
    Y = lmi(Y);
elseif isa(Y,'sdpvar')
    Y = lmi(Y);
end

% Support set+[]
if isempty(X)
    X = Y;
    return
elseif isempty(Y)   
    return
end

if ~((isa(X,'lmi')) & (isa(Y,'lmi')))
    error('Both arguments must be constraints')
end

nX = length(X.LMIid);
nY = length(Y.LMIid);
if nX==0
    X = Y;
    return
end
if nY == 0
    return;
end

xBlock = isa(X.clauses{1},'cell');
yBlock = isa(Y.clauses{1},'cell');
if yBlock && xBlock && length(X.clauses{1})>1  && length(Y.clauses{1})>1
    % Both objects are long blocks of constraints. Join these on high level
    for i = 1:length(Y.clauses)
        X.clauses{end+1} = Y.clauses{i};
    end
elseif ~xBlock && ~yBlock
    % Both have been flattened
    
    % Maybe we should build cells of cells
    if length(X.clauses)> 99
        temp = X.clauses;
        X.clauses = [];
        X.clauses{1} = temp;
        X.clauses{2} = {};
        jx = length(X.clauses);
        for i = 1:length(Y.clauses)
            X.clauses{jx}{end+1} = Y.clauses{i};
        end
    else
        for i = 1:length(Y.clauses)
            X.clauses{end+1} = Y.clauses{i};
        end
    end
else
    % This is the standard case. X has been populated with a bunch of
    % constraints (growing list) while Y is yet another constraint to be
    % added to that list.
    if ~xBlock
        % Lift to a single block
        temp = X.clauses;
        X.clauses = [];
        X.clauses{1} = temp;
    end
    % New block in X? This is where performance comes from. Handling cells
    % of cells is way faster in MATLAB, than a long cell (quadratic running
    % time, appears to be nasty copying going on when adding new elements)
    if length(X.clauses{end}) >= 100
        X.clauses{end+1} = {};
    end
    % Special case for performance
    if isa(Y.clauses{1},'cell') && length(Y.clauses)==1
        j = length(X.clauses);
        for i = 1:length(Y.clauses{1})
            X.clauses{j}{end+1} = Y.clauses{1}{i};
        end        
    else
        Y = flatten(Y);
        j = length(X.clauses);
        for i = 1:length(Y.clauses)
            X.clauses{j}{end+1} = Y.clauses{i};
        end
    end
end
aux = X.LMIid;
X.LMIid = [X.LMIid Y.LMIid];

% VERY FAST UNIQUE BECAUSE THIS IS CALLED A LOT OF TIMES....
if ~(max(aux) < min(Y.LMIid))
    i = sort(X.LMIid);
    i = i(diff([i NaN])~=0);
    if length(i)<nX+nY
        % Flatten first. This typically doesn't happen, so we accept this.
        X.clauses = [X.clauses{:}];
        [i,j] = unique(X.LMIid);
        X = subsref(X,struct('type','()','subs',{{j}}));
    end
end
