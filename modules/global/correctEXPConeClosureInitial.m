function p = correctEXPConeClosureInitial(p)
if any(p.K.e)
    % Check if we have a 0 in division. If so, just replace
    % initial with a constant and hope for the best
    top = startofEXPCone(p.K);
    % FIXME: For nonlinear cones, should work
    % with the full monomial list. Now we just crash out
    try
    for i = 1:p.K.e
        x2 = p.F_struc(top+1,:);
        if x2*[1;p.x0] == 0
            r = find(x2(2:end));
            if any(r)
                %if not any, we're doomed...
                p.x0(r(1))=1/pi;
            end
        end
        top = top + 3;
    end
    catch
    end
end