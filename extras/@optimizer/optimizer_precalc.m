function sys = optimizer_precalc(sys)

sys.model.precalc.newmonomtable = sys.model.monomtable;
sys.model.precalc.rmvmonoms = sys.model.precalc.newmonomtable(:,sys.model.parameterIndex);
sys.model.precalc.newmonomtable(:,sys.model.parameterIndex) = 0;
sys.model.precalc.Qmap = [];
% R2012b...
try
    [ii,jj,kk] = stableunique(sys.model.precalc.newmonomtable*gen_rand_hash(0,size(sys.model.precalc.newmonomtable,2),1));
    sys.model.precalc.S = sparse(kk,1:length(kk),1);
    sys.model.precalc.skipped = setdiff(1:length(kk),jj);    
    sys.model.precalc.blkOneS = blkdiag(1,sys.model.precalc.S');     
catch  
end
    
if 1%sys.nonlinear & ~sys.complicatedEvalMap
    
    % Precompute some structures
    newmonomtable = sys.model.monomtable;
    rmvmonoms = newmonomtable(:,[sys.model.parameterIndex(:);sys.model.evalParameters(:)]);
    % Linear indexation to fixed monomial terms which have to be computed
    % [ii1,jj1] = find((rmvmonoms ~= 0) & (rmvmonoms ~= 1));
    [ii1,jj1] = find( rmvmonoms < 0 | rmvmonoms > 1 | fix(rmvmonoms) ~= rmvmonoms);    
    sys.model.precalc.index1 = sub2ind(size(rmvmonoms),ii1,jj1);    
    sys.model.precalc.jj1 = jj1;    
            
    % Linear indexation to linear terms
    linterms = rmvmonoms == 1;
    if ~isempty(jj1) | any(sum(linterms,2)>1)
        [ii2,jj2] = find(linterms);
        sys.model.precalc.index2 = sub2ind(size(rmvmonoms),ii2,jj2);
        sys.model.precalc.jj2 = jj2;
        sys.model.precalc.aux = rmvmonoms*0+1;
    else
        [ii2,jj2] = find(linterms);
        sys.model.precalc.index2 = ii2;
        sys.model.precalc.jj2 = jj2;
        sys.model.precalc.aux = ones(size(rmvmonoms,1),1);
    end
    
    sys.model.newmonomtable = sys.model.monomtable;
    sys.model.rmvmonoms =  sys.model.newmonomtable(:,[sys.model.parameterIndex(:);sys.model.evalParameters(:)]);
    sys.model.newmonomtable(:,union(sys.model.parameterIndex,sys.model.evalParameters)) = 0;
   
    sys.model.removethese = find(~any(sys.model.newmonomtable,2));
    sys.model.keepingthese = find(any(sys.model.newmonomtable,2));    
end
