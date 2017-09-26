function F = subsref(F,Y)
%subsref           Overloaded indexing

switch Y(1).type
    case '()'
        thelmi = Y.subs{1};

        % Early bail to optmize code for standard case
        if isa(thelmi,'double') % &(isa(F,'lmi') %

        else
            if ~(isa(F,'lmi') & (isa(thelmi,'logical') | isa(thelmi,'double') | isa(thelmi,'char')))
                error('First argument should be an lmi object and second argument integer or string')
            end

            % If indexed using handle, flatten the structure (to simplify
            % code) and convert to index 
            if isa(thelmi,'char')
                F = flatten(F);
                thelmitemp = [];
                for i = 1:size(F.clauses,2)
                    if strcmp(F.clauses{i}.handle,thelmi)
                        thelmitemp=[thelmitemp i];
                    end
                end
                if isempty(thelmitemp)
                    em = ['LMI ''' thelmi ''' not available.'];
                    error(em)
                else
                    thelmi = thelmitemp;
                end
            end

            if islogical(thelmi)
                thelmi = double(find(thelmi));
            end

        end

        % In complicated indexing, we flatten the storing layers to
        % simplify code. 
        % For single index, we search through all buckets 
        if length(thelmi) > 1
            F = flatten(F);                   
        end
        
        if isempty(thelmi)
            F.clauses = [];
            F.LMIid = [];
            return
        end
        
        if isa(F.clauses{1},'struct') % Flat storage
            % Checks so that it exist            
            if any(thelmi<1) | any((thelmi>size(F.clauses,2)))
                em = ['LMI #' num2str(thelmi) ' not available.'];
                error(em)
            end           
            F.clauses = F.clauses(thelmi);
            F.LMIid   = F.LMIid(thelmi);           
        else        
            % Loop through all buckets
            n = cellfun(@length,F.clauses);
            cumsum_n = cumsum(n);
            bucket = min(find(thelmi <= cumsum_n));
            if isempty(bucket)                            
                found = 0;
            else
                F.LMIid   = F.LMIid(thelmi);
                if bucket > 1
                    thelmi = thelmi - cumsum_n(bucket-1);
                end               
                F.clauses = F.clauses{bucket}(thelmi);
                found = 1;
            end
%             bucket = 1;
%             top = 0;
%             found = 0;
%             while bucket <= length(F.clauses) & ~found
%                 if thelmi <= top + length(F.clauses{bucket})
%                     F.LMIid   = F.LMIid(thelmi);
%                     thelmi = thelmi - top;
%                     F.clauses = F.clauses{bucket}(thelmi);
%                     if ~isa(F.clauses,'cell')
%                         F.clauses = {F.clauses};
%                     end
%                     found = 1;
%                 else
%                     top = top  + length(F.clauses{bucket});
%                     bucket = bucket + 1;
%                 end
%             end
            if ~found
                em = ['LMI #' num2str(thelmi) ' not available.'];
                error(em)
            end                        
        end
    case '{}'
        switch length(Y)
            case 1
                if length(Y.subs)>1
                    error('Only one-dimensional indexing of set objects allowed');
                end
                if Y.subs{1}>length(F.clauses)
                    em = ['LMI #' num2str(Y.subs{1}) ' not available.'];
                    error(em)
                    error('Index exceeds set dimension.')
                end
                F = F.clauses{Y.subs{1}}.data;
            case 2
                F = F.clauses{Y(1).subs{1}}.data;
                F = subsref(F,Y(2));
            otherwise
                error('Indexing type not supported');
        end
    otherwise
        error('Indexing type not supported');
end



