function p = presolve(p)

model = p.model;

p.model.Q = p.model.Q*0;
p.model.F_struc(1:p.dimin(1),:)=[];
p.model.K.f = p.model.K.f-p.dimin(1);
redundant = sparse(zeros(p.model.K.l,1));
for i = 1:p.model.K.l
    b = p.model.F_struc(p.model.K.f+i,1);
    p.model.F_struc(p.model.K.f+i,1) = b+1;
    p.model.c = p.model.F_struc(p.model.K.f+i,2:end)';
   
     eval(['output = ' p.model.solver.call '(p.model);']);
     if output.problem == 0
         if p.model.c'*output.Primal > -b
             redundant(i) = 1
         end
     end
     
     p.model.F_struc(p.model.K.f+i,1) = b;     
end
p.model.c = c;
p.model.Q = c;