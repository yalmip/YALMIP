function p = presolve_remove_repeatedrows(p)
% Remove repeated equalities and (LP) inequalities
if p.K.l+p.K.f > 0 && p.feasible
    neq = 0;
    nin = 0;
    hash = randn(size(p.F_struc,2),1);
    hash = p.F_struc*hash;
    hash = hash(1:p.K.l+p.K.f);
    [vals, rows] = unique(hash,'stable');
    % N.B, since we do it using stable unique
    % an equality will be selected before inequality
    % and thus f(x)=0, f(x)>=0 will be presolved to f(x)=0
    remove = setdiff(1:p.K.l+p.K.f,rows);
    if ~isempty(remove)
        p.F_struc(remove,:)=[];
        neq = neq + nnz(remove > p.K.f);
        nin = nin + nnz(remove <= p.K.f);
        p.K.l = p.K.l - neq;
        p.K.f = p.K.f - nin;                
    end
    if p.K.f > 0
        % Switch the sign on the equalities, and check again
        % This can be done by simply switchg the sign on the hash value
        % Since we might have changed model since last run recompute...
        hash = randn(size(p.F_struc,2),1);
        hash = p.F_struc*hash;
        hash = hash(1:p.K.l+p.K.f);
        hash(1:p.K.f) = - hash(1:p.K.f);
        [~, rows] = unique(hash,'stable');
        remove = setdiff(1:p.K.l+p.K.f,rows);
        if ~isempty(remove)
            remove = setdiff(1:p.K.l+p.K.f,rows);
            p.F_struc(remove,:)=[];
            neq = neq + nnz(remove > p.K.f);
            nin = nin + nnz(remove <= p.K.f);
            p.K.l = p.K.l - nnz(remove > p.K.f);
            p.K.f = p.K.f - nnz(remove <= p.K.f);                        
        end
    end
    if p.options.verbose > 1 && neq+nin > 0
        disp(['* Removed ' num2str(neq+nin) ' repeated (in)equalities']);
	end
end