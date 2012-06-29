function F = subsref(F,Y)
%subsref           Overloaded indexing

% Author Johan Löfberg
% $Id: subsref.m,v 1.4 2005-02-04 10:10:27 johanl Exp $

switch Y(1).type
    case '()'
        thelmi = Y.subs{1};

        % Early bail to optmize code for standard case
        if isa(thelmi,'double') % &(isa(F,'lmi') %

        else
            if ~(isa(F,'lmi') & (isa(thelmi,'logical') | isa(thelmi,'double') | isa(thelmi,'char')))
                error('First argument should be an lmi object and second argument integer or string')
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

            if islogical(thelmi)
                thelmi = double(find(thelmi));
            end

        end

        % Checks so that it exist
        if any((thelmi<1)) | any((thelmi>size(F.clauses,2)))
            em = ['LMI #' num2str(thelmi) ' not available.'];
            error(em)
        end

        % These indicies
        j = thelmi;
        ineqs = find(j<=size(F.clauses,2));
        ineqs_j = j(ineqs);

        if isempty(ineqs_j)
            F.clauses = {};
            F.LMIid = [];
        else
            F.clauses = F.clauses(ineqs_j);
            F.LMIid   = F.LMIid(ineqs_j);
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



