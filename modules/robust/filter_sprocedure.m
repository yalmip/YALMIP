function [F_xw,F_sprocedure] = filter_sprocedure(F_xw,w,uncertaintyModel,ops)

% Select the constraints to analyse
if any(is(F_xw,'elementwise'))
    F_lp = F_xw(find(is(F_xw,'elementwise')));
    % and save the rest for later analysis
    F_xw = F_xw(find(~is(F_xw,'elementwise')));
    
    p = sdpvar(F_lp);
    keep = ones(1,length(p));
    F_sprocedure = [];
    
    % We cannot use a duality based SOS-decomposition
    opsin = ops;
    ops.verbose = 0;
    ops.sos.model=2;
    
    % All decision variables (i.e. not variables in SOS)
    Parameters = [recover(setdiff(depends(p),depends(w)))];
    for i = 1:length(p)
        d = degree(p(i),w);
        if all(d<=2) & any(d==2)
            if opsin.verbose
                disp(' - Using exact S-procedure to eliminate uncertainty');
            end
            lambda = sdpvar(1);
            if isempty(uncertaintyModel{1}.r)
                g = uncertaintyModel{1}.g;
            else
                e = (w-uncertaintyModel{1}.center);
                g = uncertaintyModel{1}.r^2-e'*e;
            end
            Parameters = [Parameters;lambda];
            F_sprocedure = [F_sprocedure, sos(p(i)-lambda*g), lambda>=0];
            %F_sprocedure = [F_sprocedure, lambda>0, compilesos(sos(p(i)-lambda*g),[],ops,s)];
            keep(i) = 0;
        elseif 0%any(d>2)
            [multiplier, lambda] = polynomial(w,2);
            e = (w-uncertaintyModel{1}.center);
            g = uncertaintyModel{1}.r^2-e'*e;
            Parameters = [Parameters;lambda];
            F_sprocedure = [F_sprocedure,lambda>=0,sos(p(i)-multiplier*g), sos(multiplier)];
            %F_sprocedure = [F_sprocedure, lambda>0, compilesos(sos(p(i)-multiplier*g)+sos(multiplier),[],ops,s)];
            keep(i) = 0;
        end
    end
    if ~isempty(F_sprocedure)
        F_sprocedure = compilesos(F_sprocedure,[],ops,Parameters);
    end
    if any(keep)
        F_xw = [F_xw,p(find(keep))>=0];
    end
else
    F_sprocedure = [];
end