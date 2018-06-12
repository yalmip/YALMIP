function pos = findhashsorted(x,searchfor)

if numel(x) < 1000
    % Short enough use built-in despite being sorted...
    pos = find(x == searchfor);
else
    % Long enough to warrant a bisection search, JITs fine
    if searchfor>x(end)
        pos = [];
        return
    elseif searchfor < x(1)
        pos = [];
        return
    elseif searchfor == x(1)
        pos = 1;
        return
    elseif searchfor == x(end)
        pos = numel(x);
        return
    end    
    low=1;
    high=numel(x);
    pos=[];
    d=numel(x);    
    while high-low > 1
        mid = floor((low+high)/2);
        val = x(mid);
        if val > searchfor
            high = mid;
        elseif val < searchfor
            low = mid;
        else
            pos = mid;
            break
        end
    end    
end