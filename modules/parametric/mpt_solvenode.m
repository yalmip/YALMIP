function model = mpt_solvenode(Matrices,lower,upper,OriginalModel,model,options)
% This is the core code. Lot of pre-processing to get rid of strange stuff
% arising from odd problems, big-M etc etc

Matrices.lb = lower;
Matrices.ub = upper;

% Remove equality constraints and trivial stuff from big-M
[equalities,redundant] = mpt_detect_fixed_rows(Matrices);
if ~isempty(equalities)
    % Constraint of the type Ex == W, i.e. lower-dimensional
    % parametric space
    if any(sum(abs(Matrices.G(equalities,:)),2)==0)
        return
    end
end
Matrices = mpt_collect_equalities(Matrices,equalities);
go_on_reducing = size(Matrices.Aeq,1)>0;
Matrices = mpt_remove_equalities(Matrices,redundant);
[Matrices,infeasible] = mpt_project_on_equality(Matrices);

% We are not interested in explicit solutions over numerically empty regions
parametric_empty = any(abs(Matrices.lb(end-Matrices.nx+1:end)-Matrices.ub(end-Matrices.nx+1:end)) < 1e-6);

% Were the equality constraints found to be infeasible?
if infeasible | parametric_empty
    return
end

% For some models with a lot of big-M stuff etc, the amount of implicit
% equalities are typically large, making the LP solvers unstable if they
% are not removed. To avoid problems, we iteratively detect fixed variables
% and strenghten the bounds.
fixed = find(Matrices.lb == Matrices.ub);
infeasible = 0;
while 0%~infeasible & options.mp.presolve
    [Matrices,infeasible] = mpt_derive_bounds(Matrices,options);
    if isequal(find(Matrices.lb == Matrices.ub),fixed)
        break
    end
    fixed = find(Matrices.lb == Matrices.ub);
end

% We are not interested in explicit solutions over numerically empty regions
parametric_empty = any(abs(Matrices.lb(end-Matrices.nx+1:end)-Matrices.ub(end-Matrices.nx+1:end)) < 1e-6);

if ~infeasible & ~parametric_empty

    while go_on_reducing & ~infeasible & options.mp.presolve
        [equalities,redundant] = mpt_detect_fixed_rows(Matrices);
        if ~isempty(equalities)
            % Constraint of the type Ex == W, i.e. lower-dimensional
            % parametric space
            if any(sum(abs(Matrices.G(equalities,:)),2)==0)
                return
            end
        end
        Matrices = mpt_collect_equalities(Matrices,equalities);
        go_on_reducing = size(Matrices.Aeq,1)>0;
        Matrices = mpt_remove_equalities(Matrices,redundant);
        [Matrices,infeasible] = mpt_project_on_equality(Matrices);
        M = Matrices;
        if  go_on_reducing & ~infeasible
            [Matrices,infeasible] = mpt_derive_bounds(Matrices,options);
        end
        if infeasible
            % Numerical problems most likely because this infeasibility
            % should have been caught above. We have only cleaned the model
            Matrices = M;
        end
    end

    if ~infeasible
        if Matrices.qp
            if isequal(options.solver,'mplcp')
                [Pn,Fi,Gi,ac,Pfinal,details] = mpt_mpqp_mplcp(Matrices,options);
            else
                e = eig(full(Matrices.H));
                if min(e) == 0
                    disp('Lack of strict convexity may lead to troubles in the mpQP solver')
                elseif min(e) < -1e-8
                    disp('Problem is not positive semidefinite! Your mpQP solution may be completely wrong')
                elseif min(e) < 1e-5
                    disp('QP is close to being (negative) semidefinite, may lead to troubles in mpQP solver')
                end
                %Matrices.H = Matrices.H + eye(length(Matrices.H))*1e-4;
                [Pn,Fi,Gi,ac,Pfinal,details] = mpt_mpqp(Matrices,options.mpt);
            end
        else
            if isequal(options.solver,'mplcp')
                [Pn,Fi,Gi,ac,Pfinal,details] = mpt_mpqp_mplcp(Matrices,options);
            else               
                [Pn,Fi,Gi,ac,Pfinal,details] = mpt_mplp(Matrices,options.mpt);
            end
        end
        [Fi,Gi,details] = mpt_project_back_equality(Matrices,Fi,Gi,details,OriginalModel);
        [Fi,Gi] = mpt_select_rows(Fi,Gi,Matrices.requested_variables);
        [Fi,Gi] = mpt_clean_optmizer(Fi,Gi);
        model = mpt_appendmodel(model,Pfinal,Pn,Fi,Gi,details);
       % model = mpt_reduceOverlaps_orderfaces(model);if ~isa(model,'cell');model = {model};end
    end
else

end

function [Pn,Fi,Gi,ac,Pfinal,details] = mpt_mpqp_mplcp(Matrices,options)

if Matrices.qp
    lcpData = lcp_mpqp(Matrices);
    BB = mplcp(lcpData)
    [Pn,Fi,Gi] = soln_to_mpt(lcpData,BB);
else
    lcpData = lcp_mplp(Matrices);
    BB = mplcp(lcpData)
    [Pn,Fi,Gi] = soln_to_mpt(lcpData,BB);
end
Pfinal = union(Pn);
if Matrices.qp
    for i=1:length(Fi)
        details.Ai{i} = 0.5*Fi{i}'*Matrices.H*Fi{i} + 0.5*(Matrices.F*Fi{i}+Fi{i}'*Matrices.F') + Matrices.Y;
        details.Bi{i} = Matrices.Cf*Fi{i}+Gi{i}'*Matrices.F' + Gi{i}'*Matrices.H*Fi{i} + Matrices.Cx;
        details.Ci{i} = Matrices.Cf*Gi{i}+0.5*Gi{i}'*Matrices.H*Gi{i} + Matrices.Cc;
    end
else
    for i=1:length(Fi)
        details.Ai{i} = [];
        details.Bi{i} = Matrices.H*Fi{i};
        details.Ci{i} = Matrices.H*Gi{i};
    end
end
ac  = [];
