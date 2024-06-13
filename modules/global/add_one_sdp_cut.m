function  [X,p_lp,infeasibility,asave,bsave,failure] = add_one_sdp_cut(p,p_lp,x,i,p_original)
% Common function for CUTSDP and BNB to generate a cut at infeasible point
newcuts = 0;
newF = [];
n = p.K.s(i);
if numel(x)/length(x) < .1
    X = p.semidefinite{i}.F_struc*sparse([1;x]);
else
    X = p.semidefinite{i}.F_struc*[1;x];
end
X = reshape(X,n,n);X = (X+X')/2;
asave = [];
bsave = [];
% First check if it happens to be psd. Then we are done. Quicker
% than computing all eigenvalues
% This also acts as a slight safe-guard in case the sparse eigs
% fails to prove that the smallest eigenvalue is non-negative
%[R,indefinite] = chol(X+eye(length(X))*1e-12);
%if indefinite

% User is trying to solve by only generating no-good cuts
permutation = [];
failure = 0;
if p.options.cutsdp.cutlimit == 0
     [d,v] = eig(full(X));
     infeasibility = v(1,1);   
    return
end

% For not too large problems, we simply go with a dense
% eigenvalue/vector computation
if  n <= p_lp.options.cutsdp.switchtosparse
    [d,v] = eig(full(X));
    failure = 0;
else
    % Try to perform a block-diagonalization of the current solution,
    % and compute eigenvalue/vectors for each block.
    % Sparse eigenvalues can easily fail so we catch info about this
    [d,v,permutation,failure] = dmpermblockeig(X,p_lp.options.cutsdp.switchtosparse);
end
d(abs(d)<1e-12)=0;
% This indicator is used by CUTSDP to judge feasibility
infeasibility = min(diag(v));
% However, we might be interested in cuts also for feasible case
% when calling this to generate cuts in BNB
if infeasibility<=1e-6
    [ii,jj] = sort(diag(v));
    
    if ~isempty(permutation)
        [~,inversepermutation] = ismember(1:length(permutation),permutation);      
    end
    
    for m = jj(1:min(length(jj),p.options.cutsdp.cutlimit))'
        if v(m,m)<=1e-5                   
            try
                if ~isempty(permutation)
                    dhere = d(inversepermutation,m);
                else
                    dhere = d(:,m);
                end
                dd = dhere*dhere';dd = dd(:);
                bA = dd'*p.semidefinite{i}.F_struc;
                if numel(bA)/nnz(bA) < .1
                    bA = sparse(bA);
                end
            end
            b = bA(:,1);
            A = -bA(:,2:end); 
            if isempty(p_lp.F_struc) || ~any(sum(abs(p_lp.F_struc-repmat([b -A],size(p_lp.F_struc,1),1)),2)<= 1e-12)
                newF = real([newF;[b -A]]);
                newcuts = newcuts + 1;
                if isempty(asave)
                    A(abs(A)<1e-12)=0;
                    b(abs(b)<1e-12)=0;
                    asave = -A(:);
                    bsave = b;
                end
            end
        end
    end
end
newF(abs(newF)<1e-12) = 0;
keep=find(any(newF(:,2:end),2));
newF = newF(keep,:);
if size(newF,1)>0
    newF(:,1) = newF(:,1) + 0*0.02*abs(p_lp.options.cutsdp.feastol);
    p_lp.F_struc = [p_lp.F_struc(1:p_lp.K.f,:);p_lp.F_struc(1+p_lp.K.f:end,:);newF];
    p_lp.K.l = p_lp.K.l + size(newF,1);
end

