function pos = findhashsorted(T,t)

% This code does not exolit the sorted hash. It is only exploited in the
% assocoated c-version
if isempty(T)
    pos = []; % Random warnings on 6.1
else
    pos = find(T==t);
end
    