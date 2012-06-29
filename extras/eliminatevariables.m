function [model,keptvariables,infeasible] = eliminatevariable(model,varindex,value)

keptvariables = 1:length(model.c);

model.F_struc(1:length(varindex),:) = [];
model.K.f = model.K.f - length(varindex);

newmonomtable = model.monomtable;
rmvmonoms = newmonomtable(:,varindex);
newmonomtable(:,varindex) = 0;

monomvalue = prod(repmat(value(:)',size(rmvmonoms,1),1).^rmvmonoms,2);

removethese = find(~any(newmonomtable,2));
keepingthese = find(any(newmonomtable,2));

value = monomvalue(removethese);
monomgain = monomvalue;monomgain(removethese) = [];

if ~isempty(model.F_struc)
    model.F_struc(:,1) = model.F_struc(:,1)+model.F_struc(:,1+removethese)*value;
    model.F_struc(:,1+removethese) = [];
    model.F_struc = model.F_struc*diag([1;monomgain]);
end
infeasible = 0;
if model.K.f > 0
    candidates = find(sum(abs(model.F_struc(1:model.K.f,2:end)),2) == 0);
    if ~isempty(candidates)
        infeasibles = find(model.F_struc(candidates,1)~=0);
        if ~isempty(infeasibles)
            infeasible = 1;
            return
        else
            model.F_struc(candidates,:) = [];
            model.K.f = model.K.f - length(candidates);
        end
    end
end
if model.K.l > 0
    candidates = find(sum(abs(model.F_struc(model.K.f + (1:model.K.l),2:end)),2) == 0);
    if ~isempty(candidates)
        infeasibles = find(model.F_struc(model.K.f + candidates,1)<0);
        if ~isempty(infeasibles)
            infeasible = 1;
            return
        else
            model.F_struc(model.K.f + candidates,:) = [];
            model.K.l = model.K.l - length(candidates);
        end
    end
end

model.f = model.f + model.c(removethese)'*value;

model.c(removethese)=[];
model.c = model.c + 2*model.Q(keepingthese,removethese)*value;
model.Q(removethese,:) = [];
model.Q(:,removethese) = [];

model.c = model.c.*monomgain;
keptvariables(removethese) = [];

model.lb(removethese)=[];
model.ub(removethese)=[];

newmonomtable(:,removethese) = [];
newmonomtable(removethese,:) = [];

skipped = [];
for i  = 1:size(newmonomtable,1)
    this = newmonomtable(i,:);
    j = findrows(newmonomtable(i+1:1:end,:),this);
    if ~isempty(j)
        j = j + i;
        j = setdiff(j,skipped);
        if ~isempty(j)
            model.c(i) = model.c(i) + sum(model.c(j));
            model.F_struc(:,i+1) = model.F_struc(:,i+1) + sum(model.F_struc(:,j+1),2);
            skipped = [skipped j(:)'];
        end
    end
end
model.c(skipped) = [];
model.lb(skipped) = [];
model.ub(skipped) = [];
model.F_struc(:,1+skipped) = [];
newmonomtable(skipped,:) = [];
newmonomtable(:,skipped) = [];
model.Q(:,skipped)=[];
model.Q(skipped,:)=[];
keptvariables(skipped) = [];

model.monomtable = newmonomtable;

model.variabletype = zeros(size(model.monomtable,1),1)';
nonlinear = ~(sum(model.monomtable,2)==1 & sum(model.monomtable~=0,2)==1);
if ~isempty(nonlinear)
    model.variabletype(nonlinear) = 3;
    quadratic = sum(model.monomtable,2)==2;
    model.variabletype(quadratic) = 2;
    bilinear = max(model.monomtable,[],2)<=1;
    model.variabletype(bilinear & quadratic) = 1;
    sigmonial = any(0>model.monomtable,2) | any(model.monomtable-fix(model.monomtable),2);
    model.variabletype(sigmonial) = 4;
end

model.x0 = zeros(length(model.c),1);
%model.x0 = zeros(length(find(model.variabletype)),1);