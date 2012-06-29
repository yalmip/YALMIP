function [oldvariable,sqrList] = expandforindex(sqrList,candidates,nonconvindex);
%EXPANDFORINDEX Internal function for organizing nonlinear variables

% Author Johan Löfberg
% $Id: expandforindex.m,v 1.4 2004-08-06 09:02:52 johanl Exp $


% Expand the incoming variables

if 1
z = sqrList(:,1);

x1_list = find(z==candidates(2));
if candidates(2)==candidates(1)
    x2_list = x1_list;
else
    x2_list = find(z==candidates(1));
end

temp = [];
if isempty(x1_list)
    temp = candidates(2);
else
    temp = sqrList(x1_list,2:end);
end
if isempty(x2_list)
    temp = [temp candidates(1)];
else
    temp = [temp sqrList(x2_list,2:end)];
end
temp = sort(temp(find(any(temp,1))));

oldsqrList = sqrList;
if length(temp)>size(sqrList,2)-1
    % This one cannot already exist
    index = [];
    
    %sqrList = [sqrList zeros(size(sqrList,1),1+size(sqrList,2)-length(temp));nonconvindex temp(end:-1:1)];   
    sqrList = [sqrList zeros(size(sqrList,1),1+length(temp)-size(sqrList,2));nonconvindex temp(end:-1:1)];   
else
%    sqrList = [sqrList;nonconvindex temp(end:-1:1) zeros(1,size(sqrList,2)-length(temp)-1)];
%    searchfor = sqrList(end,2:end);

    searchfor = [temp(end:-1:1) zeros(1,size(sqrList,2)-length(temp)-1)];
    searchin = sqrList(:,2:end);
    
    % Simple Hash
    key = sum(searchfor);
    tbl = sum(searchin,2);
    
    possible = find(tbl==key);
    
    index = findrows(searchin(possible,:),searchfor);
    if isempty(index)
         sqrList = [sqrList;nonconvindex temp(end:-1:1) zeros(1,size(sqrList,2)-length(temp)-1)];
    else
         index = possible(index);
    end
end

%searchfor = sqrList(end,2:end);
%index = findrows(sqrList(1:end-1,2:end),searchfor);

if length(index)>1
    index
end

if ~isempty(index)
    oldvariable = sqrList(index,1);
    sqrList = oldsqrList;
else
    oldvariable = [];
end



else


% Expand the list to begin with
if isempty(sqrList)
    sqrList = [nonconvindex candidates];
else
    sqrList = [sqrList; nonconvindex candidates(2) candidates(1) zeros(1,size(sqrList,2)-3)];        
end
bottom = size(sqrList,1); 

%%FIXz = sqrList(1:bottom,1);
z = sqrList(:,1);

x1_list = find(z==candidates(2));
if candidates(2)==candidates(1)
    x2_list = x1_list;
else
    x2_list = find(z==candidates(1));
end

temp = [];
if isempty(x1_list)
    temp = sqrList(bottom,2);
else
    temp = sqrList(x1_list,2:end);
end
if isempty(x2_list)
    temp = [temp sqrList(bottom,3)];
else
    temp = [temp sqrList(x2_list,2:end)];
end
temp = sort(temp(find(any(temp,1))));
sqrList(bottom,2:2+length(temp)-1)=temp(end:-1:1);%fliplr(sort(temp));

bottom = bottom+1;

searchfor = sqrList(end,2:end);
index = findrows(sqrList(1:end-1,2:end),searchfor);

if ~isempty(index)
    oldvariable = sqrList(index,1);
    sqrList = sqrList(1:end-1,:);
else
    oldvariable = [];
end

end
