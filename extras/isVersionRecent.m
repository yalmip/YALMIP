function isRecent = isVersionRecent(version1, version2)
% Author: Chat-GPT

% Split the version strings into arrays of integers
version1Arr = str2double(strsplit(version1, '.'));
version2Arr = str2double(strsplit(version2, '.'));

% Add zero components to version arrays if they are not in the format x.y.z
if numel(version1Arr) < 3
    version1Arr(end+1:3) = 0;
end
if numel(version2Arr) < 3
    version2Arr(end+1:3) = 0;
end

% Compare the version numbers
if version1Arr(1) > version2Arr(1)
    isRecent = true;
elseif version1Arr(1) < version2Arr(1)
    isRecent = false;
elseif version1Arr(2) > version2Arr(2)
    isRecent = true;
elseif version1Arr(2) < version2Arr(2)
    isRecent = false;
elseif version1Arr(3) >= version2Arr(3)
    isRecent = true;
else % Version 1 is older than Version 2
    isRecent = false;
end
end
