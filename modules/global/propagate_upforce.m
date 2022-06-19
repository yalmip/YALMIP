function [p,p_feasible] = propagate_upforce(p)
p_feasible = 1;

for i = 1:length(p.upForce)
    forcing = p.upForce{i}.forcing;
    forced = p.upForce{i}.forced;
    if p.lb(forcing)==1
        [p.lb(forced) p.ub(forced)]
%         if any(p.lb(forced)==1)
%             error
%         else
%          %   p.ub(forced)=0;
%         end
    end
end
    