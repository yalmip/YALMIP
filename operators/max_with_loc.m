function [y,loc] = max_with_loc(varargin)

% Called from max to handle case when user wants location index

X = varargin{1};

% Trivial case
if nargin == 3 && varargin{3}==2 && size(varargin{1},2)==1
    y = varargin{1};
    loc = ones(length(y),1);
    return
end
if nargin == 3 && varargin{3}==1 && size(varargin{1},1)==1
    y = varargin{1};
    loc = ones(1,length(y));
    return
end

% first simple case, simply vector argument
if min(size(X))==1
    [temp1,temp2] = sort(X);
    y = temp1(end);
    loc = temp2(end);
    return
end

% case 2, max(X)
if nargin == 1 && min(size(X))>1
    y = [];
    loc = [];
    for i = 1:size(X,2)
         [temp1,temp2] = sort(X(:,i));
         y = [y temp1(end)];
         loc = [loc temp2(end)];
    end
    return
end

% Specified direction
if nargin == 3
    if varargin{3}==1
        y = [];
        loc = [];
        for i = 1:size(X,2)
            [temp1,temp2] = sort(X(:,i));
            y = [y temp1(end)];
            loc = [loc temp2(end)];
        end
    else
        [y,loc] = max(X');
        y = y';
        loc = loc';
    end
    return
end

