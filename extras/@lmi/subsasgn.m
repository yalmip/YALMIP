function F0 = subsasgn(F,Y,Z)
%subasgn           Overloaded indexing

% Author Johan Löfberg 
% $Id: subsasgn.m,v 1.4 2005-04-29 16:28:03 joloef Exp $   

% TODO : Implement lmi/subsasgn...
%error('set/subsasgn not fully implemented yet. Contact author if you need this functionality')
switch Y(1).type
    case '()'
        thelmi = Y.subs{1};
        
        if ~isa(F,'lmi')
            error('First argument should be a SET object')
        end
        
        if ~(isa(thelmi,'double') | isa(thelmi,'char'))
            error('Index should be an integer integer or tag-string')
        end
        
        if ~(isa(Z,'lmi') | isempty(Z))
            error('Right-hand side should be a SET object or empty');
        end
        
        % If indexed using handle, convert to index
        if isa(thelmi,'char')
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
        
        % Convert to linear indecicies
        if islogical(thelmi)
           thelmi = double(find(thelmi));
        end
        
        % Checks so that it exist
        if any((thelmi<1)) | any((thelmi>size(F.clauses,2)))
            em = ['LMI #' num2str(thelmi(find(thelmi>size(F.clauses,2)))) ' not available.'];
            error(em)
        end
        
        F0 = F;
        
        % Specialized for deleting constraints
        if isempty(Z)
            F0.clauses = {F0.clauses{setdiff(1:size(F.clauses,2),thelmi)}};
            F0.LMIid = F0.LMIid(setdiff(1:size(F.clauses,2),thelmi));
        else
            if length(thelmi) < length(Z.clauses)
                error('???  In an assignment  A(I) = B, the number of elements in B and I must be the same.');
            end
            for i = 1:length(thelmi)
                F0.clauses{thelmi(i)} = Z.clauses{i};
                F0.LMIid(thelmi(i)) = Z.LMIid(i);
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
                F0 = F.clauses{Y.subs{1}}.data;                
            case 2
                F0 = F.clauses{Y(1).subs{1}}.data;
                F0 = subsref(F0,Y(2));
            otherwise
                error('Indexing type not supported');   
        end
    otherwise
        error('Indexing type not supported');        
end



