function [p_sdpknapsackcuts,p_rawcuts,p_lp_contcut] = create_sdpknapsack_cuts(p,x_trial,poriginal,upper)

p_sdpknapsackcuts = emptyNumericalModel;
p_rawcuts = emptyNumericalModel;
p_sdpcuts = emptyNumericalModel;
x_trial(x_trial < p.lb) = p.lb(x_trial < p.lb);
x_trial(x_trial > p.ub) = p.ub(x_trial > p.ub);

for sdpcount = 1:length(p.K.s)
    ptemp = emptyNumericalModel;
    ptemp.options = p.options;
    [X,p_lp_contcut] = add_one_sdp_cut(p,ptemp,x_trial,sdpcount,p);
    if ~isempty(p_lp_contcut.F_struc)
        for j = 1:p_lp_contcut.K.l
            row = p_lp_contcut.F_struc(j,:);
            % Set continuous to best possible
            row = fix_continuous_in_row_at_best_possible(row,poriginal,upper);                        
            % Eigenvalues are easily somewhat shaky numerically so
            % relax slightly
            row(1) = row(1) + 1e-8; 
            % Save for furter use later
            p_rawcuts = addInequality(p_rawcuts,row);
            
            % Use general cover separator
            cut = knapsack_create_cover_cut(-row(2:end),row(1),x_trial,'gu',p.gubs);
            cut2 = knapsack_cheapcut(row);            
            ptemp.F_struc = [ptemp.F_struc;cut;cut2];
            ptemp.K.l = ptemp.K.l + size(cut,1)+size(cut2,1);            
        end
        if ptemp.K.l>0
            %violations = ptemp.F_struc*[ones(1,size(allRelaxedSolutions,2));allRelaxedSolutions];
            violations = ptemp.F_struc*[1;x_trial];
            t = find(violations(:,end) < 0);
            [nv,loc] = max(sum(violations < -1e-1,2));
            if nv > 0 && violations(loc,end) < -1e-1
                % covercandidates = [covercandidates;p_lp.F_struc(loc,:)];
                % covers = [covers;p_lp.F_struc(loc,:)];
            else
            end
            p_sdpknapsackcuts = addInequality(p_sdpknapsackcuts,ptemp.F_struc(t,:));
            %covers = [covers;ptemp.F_struc(t,:)];
        end
    end
end