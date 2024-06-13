function [sol,m,Q,residuals,everything] = solvesos_find_blocks(F,obj,options,params,candidateMonomials)

tol = options.sos.numblkdg;
if tol > 1e-2
    disp(' ');
    disp('-> Are you sure you meant to have a tolerance in numblk that big!')
    disp('-> The options numblkdiag controls the tolerance, it is not a 0/1 switch.')
    disp(' ');
end
options.sos.numblkdg = 0;
[sol,m,Q,residuals,everything] = solvesos(F,obj,options,params,candidateMonomials);

% Save old structure to find out when we have stalled
for i = 1:length(Q)
    oldlengths{i} = length(Q{i});
end

go_on = (sol.problem == 0 | sol.problem == 4);
iterations = 0;
while go_on && iterations < options.sos.numblkiterlimit
    iterations = iterations + 1;
    for sosfun = 1:length(Q)
        Qtemp = Q{sosfun};
        keep = diag(Qtemp)>tol;
        Qtemp(:,find(~keep)) = [];
        Qtemp(find(~keep),:) = [];

        m{sosfun} = m{sosfun}(find(keep));

        Qtemp(abs(Qtemp) < tol) = 0;
        [v1,dummy1,r1,dummy3]=dmperm(Qtemp+eye(length(Qtemp)));
        lengths{sosfun} = [];
        n{sosfun} = {};
        for blocks = 1:length(r1)-1
            i1 = r1(blocks);
            i2 = r1(blocks+1)-1;
            if i2>i1
                n{sosfun}{blocks} = m{sosfun}(v1(i1:i2));
            else
                n{sosfun}{blocks} = m{sosfun}(v1(i1));
            end
            lengths{sosfun} =  [lengths{sosfun}  length(n{sosfun}{blocks})];
        end
        lengths{sosfun} = sort(lengths{sosfun});
    end
    go_on = ~isequal(lengths,oldlengths);
    oldlengths = lengths;
    if go_on
        [sol,m,Q,residuals,everything] = solvesos(F,obj,options,params,n);
        go_on = go_on & (sol.problem == 0 | sol.problem == 4);
        if sol.problem == 1
            disp('-> Feasibility was lost during the numerical block-diagonalization.')
            disp('-> The setting sos.numblkdiag is probably too big')
        end
    end
end