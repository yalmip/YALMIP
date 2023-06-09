function plotNodeModel(p,color,shade)

P = [p.lpcuts;
    p.F_struc(p.K.f+1:p.K.f+p.K.l,:)]*[1;p.plotter.x]>=0;
P = [P, p.lb<=p.plotter.x<=p.ub];
if any(p.c)
    % Could be feasibility problem
    P = [P, p.c'*p.plotter.x<=p.upper];
end
if any(p.K.q)
    for i = 1:length(p.K.q)
        top = startofSOCPCone(p.K);
        e = p.F_struc(top:top+p.K.q(i)-1,:)*[1;p.plotter.x];
        P = P + cone(e);
        top = top + p.K.q(i);
    end
end
if nargin < 2
    color = 'b';
end
if nargin >= 3
    p.plotter.ops.plot.shade = shade;
end
hold on
plot(P,p.plotter.x(1:2),color,[],p.plotter.ops)