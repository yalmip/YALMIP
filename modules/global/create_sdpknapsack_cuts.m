function [p_sdpknapsackcuts,p_binarycut,p_lp_contcut] = create_sdpknapsack_cuts(p,x_trial,poriginal,upper)

p_sdpknapsackcuts = emptyNumericalModel;
p_binarycut = emptyNumericalModel;
p_lp_contcut = emptyNumericalModel;
x_trial(x_trial < p.lb) = p.lb(x_trial < p.lb);
x_trial(x_trial > p.ub) = p.ub(x_trial > p.ub);

for sdpcount = 1:length(p.K.s)
    ptemp = emptyNumericalModel;
    ptemp.options = p.options;   
    [X,cone_cuts] = add_one_sdp_cut(p,ptemp,x_trial,sdpcount,p);
    if ~isempty(cone_cuts.F_struc)
        p_lp_contcut = addInequality(p_lp_contcut,cone_cuts.F_struc);
        for j = 1:cone_cuts.K.l
            row = cone_cuts.F_struc(j,:);
            % Set continuous to best possible
            row = fix_continuous_in_row_at_best_possible(row,poriginal,upper);                        
            % Eigenvalues are easily somewhat shaky numerically so
            % relax slightly
            row(1) = row(1) + 1e-8; 
            used = find(row(2:end));
            if isempty(used) || ~all(p.isbinary(used))
                continue
            end
            % Strengthen coefficients
             if all(row(2:end) >= 0) && row(1) < 0
                 row(2:end) = min(row(2:end),-row(1));
             end
             % Trivial fix rounding 
             if all(row(2:end)==round(row(2:end)))
                 row(1) = floor(row(1));
             end
            % Save for furter use later
            p_binarycut = addInequality(p_binarycut,row);
            
            % Use general cover separator
            cut = knapsack_create_cover_cut(-row(2:end),row(1),x_trial,'gu');
            cut2 = knapsack_cheapcut(row);            
            ptemp.F_struc = [ptemp.F_struc;cut;cut2];
            ptemp.K.l = ptemp.K.l + size(cut,1)+size(cut2,1);            
        end
        if ptemp.K.l>0
            p_sdpknapsackcuts = addInequality(p_sdpknapsackcuts,ptemp.F_struc);           
        end
    end
end
