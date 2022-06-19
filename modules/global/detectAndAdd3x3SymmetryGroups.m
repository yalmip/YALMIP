function p = detectAndAdd3x3SymmetryGroups(p)

if ~any(p.K.f)
    % Put in format for the detector
    top = startofSDPCone(p.K);
    p.F_struc = p.F_struc';
    for i = 1:length(p.K.s)
        p.semidefinite{i}.F_struc = p.F_struc(:,top:top+p.K.s(i)^2-1)';
        top = top + p.K.s(i)^2;
    end
    p.F_struc = p.F_struc';
    p.sdpsymmetry = [];
    % Find
    p = detect3x3SymmetryGroups(p);
    % Create
    pp = p;pp.F_struc = [];pp.K.l = 0;pp.K.f = 0;pp.K.s = 0;
    [p,pp] = addSymmetryCuts(p,pp);
    % Add
    p = addInequality(p,pp.F_struc);
    p.semidefinite=[];
end