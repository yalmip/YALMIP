function p_lp = add3x3sdpsymmetrycut(p,p_lp,x)

for j = 1:length(p.sdpsymmetry)
    excludes = [];
    n = sqrt(size(p.sdpsymmetry{j}.dataBlock,1));
    for i = 1:length(p.sdpsymmetry{j}.variables)
        if min(eig(reshape(p.sdpsymmetry{j}.dataBlock*[1;x(p.sdpsymmetry{j}.variables{i})],n,n))) < -abs(p_lp.options.bnb.feastol)
            excludes = [excludes x(p.sdpsymmetry{j}.variables{i})];
        end
    end
    if ~isempty(excludes)
        newF = [];
        infeasible_combinations = unique(excludes','rows')';
        for k = 1:size(infeasible_combinations,2)
            % Local cut for reduced set
            [b,atemp] = exclusionCut(infeasible_combinations(:,k),-1);
            % Add that cut for every variable groups
            for s = 1:length(p.sdpsymmetry{j}.variables)
                a = spalloc(1,length(p_lp.c),1);
                a(p.sdpsymmetry{j}.variables{s}) = atemp;
                if all(sum(abs(p_lp.F_struc - [b a]),2)<=1e-12)
                    newF = [newF;b a];
                end
            end
        end
        p_lp = addInequality(p_lp,newF);
    end
end