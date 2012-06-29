function [Fi,Gi] = mpt_clean_optmizer(Fi,Gi);
if length(Fi)>0
    for i = 1:length(Fi)
        Fi{i} = round(1e10*Fi{i})/1e10;
        Gi{i} = round(1e10*Gi{i})/1e10;
    end
end