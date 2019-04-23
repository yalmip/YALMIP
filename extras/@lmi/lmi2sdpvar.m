function X = lmi2sdpvar(F,thelmi)

if isa(F.clauses{1},'struct') % Flat storage
    % Checks so that it exist
    if any(thelmi<1) | any((thelmi>size(F.clauses,2)))
        em = ['LMI #' num2str(thelmi) ' not available.'];
        error(em)
    end
    X = F.clauses{thelmi}.data;    
else
    % Loop through all buckets
    n = cellfun(@length,F.clauses);
    cumsum_n = cumsum(n);
    bucket = min(find(thelmi <= cumsum_n));
    if isempty(bucket)
        found = 0;
    else      
        if bucket > 1
            thelmi = thelmi - cumsum_n(bucket-1);
        end
        X = F.clauses{bucket}{thelmi}.data;
        found = 1;
    end
    if ~found
        em = ['LMI #' num2str(thelmi) ' not available.'];
        error(em)
    end
end




