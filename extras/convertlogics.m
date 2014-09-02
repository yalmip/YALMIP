function [F,changed] = convertlogics(F)
%CONVERTLOGICS Internal function to convert logic constraints to mixed integer constraints

changed = 0;
if length(F)>0
    extvariables = yalmip('logicextvariables');
    if ~isempty(extvariables)       
        for i = 1:length(F)
            if is(F(i),'elementwise')
                Fi = sdpvar(F(i));
                Fv =getvariables(Fi);
                if length(Fv)==1
                    xb = getbase(Fi);
                    if isequal(xb,[0 1])
                        if ismember(Fv,extvariables)
                            F(i) = (Fi >= 1);
                        end
                    end
                end
            end
        end
    end
end
