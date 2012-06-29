function cut = disjunctivecut(varargin)
% 
% 
% x = sdpvar(2,1);
% A1 = randn(8,2);
% b1 = rand(8,1)*2-A1*[3;3];
% A2 = randn(8,2);
% b2 = rand(8,1)*2-A2*[-3;3];
% A3 = randn(8,2);
% b3 = rand(8,1)*2-A3*[3;-3];
% C1 = A1*x < b1;
% C2 = A2*x < b2;
% C3 = A3*x < b3;
% cut = disjunctivecut(C1,C2,C3,[-8;0]);cut = disjunctivecut(C1,C2,C3,[-8;0]);
% plot([cut,-10<x<10]);plot(hull(C1,C2,C3),x);plot(C1,[],'y');plot(C2,[],'y');plot(C3,[],'y')
xstar = varargin{end};
x = [];
beta = sdpvar(1);
alpha = sdpvar(length(xstar),1);
for i = 1:nargin-1  
    [Imodel,Iax1,Iax2,p{i}] = export(varargin{i},[],sdpsettings,[],[],0);
    neq(i) = p{i}.K.f;
    b{i} = -p{i}.F_struc(:,1);
    A{i} = p{i}.F_struc(:,2:end);
    mu{i} = sdpvar(length(b{i}),1);
end

Objective  = alpha'*xstar-beta;
Constraints = [-1<alpha<1];

summu = 0;
for i = 1:nargin-1
    summu = summu + sum(mu{i});
    Constraints = [Constraints,alpha' == mu{i}'*A{i}];
    Constraints = [Constraints,beta  <= mu{i}'*b{i}];
    Constraints = [Constraints,mu{i}(p{i}.K.f+1:end)>0];
end

solvesdp(Constraints,Objective,sdpsettings('verbose',0));

x = recover(p{1}.used_variables);
cut = set(-double(beta)+double(alpha)'*x > 0);


return

%x = recover(sdpvar(C1));

%[Imodel,Iax1,Iax2,p1] = export(C1,[],sdpsettings,[],[],0);
%[Omodel,Oax1,Oax2,p2] = export(C2,[],sdpsettings,[],[],0);

%neq1 = p1.K.f;
%neq2 = p2.K.f;

%e1 = p1.F_struc*[1;x];
%e2 = p2.F_struc*[1;x];
%Model1 = [e2(1+p1.K.f:end)>=0];
%Model2 = [e2(1+p2.K.f:end)>=0];
if 0
Ab1 = getbase(sdpvar(Model1));
Ab2 = getbase(sdpvar(Model2));
b1 = -Ab1(:,1);
A1 =  Ab1(:,2:end);
b2 = -Ab2(:,1);
A2 = Ab2(:,2:end);

b1 = -p1.F_struc(:,1);
A1 = p1.F_struc(:,2:end);
b2 = -p2.F_struc(:,1);
A2 = p2.F_struc(:,2:end);
end

alpha = sdpvar(length(xstar),1);
beta = sdpvar(1);
mu1 = sdpvar(length(b1),1);
mu2 = sdpvar(length(b2),1);

Objective  = alpha'*xstar-beta;
Constraint = [alpha' == mu1'*A1,alpha' == mu2'*A2, beta <= mu1'*b1, beta <= mu2'*b2,mu1(neq1+1:end)>0,mu2(neq2+1:end)>0];
Constraint = [Constraint,-1<alpha<1,sum(mu1)+sum(mu2)<10];

solvesdp(Constraint,Objective,sdpsettings('verbose',0));

cut = set(-double(beta)+double(alpha)'*x > 0);
