function test_chance_2

yalmip('clear')
sdpvar a b t
sdpvar t

Model = [probability(a >= t) >= 0.5,
         probability(b >= t) >= 0.5, 
         uncertain(a,'normal',[0],eye(1)),
         uncertain(b,'normal',[0],eye(1))];
solvesdp(derandomize(Model),-t)
mbg_asserttolequal(double(t), 0, 1e-5);


yalmip('clear')
sdpvar a b t m
sdpvar t

Model = [probability(a >= t) >= 0.5,
         probability(b >= t) >= 0.5, 
         uncertain(a,'normal',[0],eye(1)),
         uncertain(b,'normal',m,eye(1)),
         uncertain(m,'normal',0,eye(1)), ];
P=optimizer(derandomize(Model),-t,sdpsettings('solver','mosek'),m,t)
Q=sample(P,20)

sdpvar w wmean
optimize(derandomize([uncertain(w,'normal',wmean,1), probability(w >= 0) >= 0.9]),wmean,sdpsettings('solver','cplex','debug',1))
mbg_asserttolequal(double(wmean), 1.2816, 1e-3);

clear all;yalmip('clear');
sdpvar w(2,1) wmean(2,1) gamma
uu = [];
for gamma = .1:.05:.9
    Model = derandomize([uncertain(w,'normal',wmean,1), probability(sum(w) >= 0) >= gamma]);
    optimize(Model,wmean'*wmean-gamma,sdpsettings('solver','','debug',1));
    uu = [uu value(wmean'*wmean-gamma)];
end

Model = [uncertain(w,'normal',wmean,1), probability(sum(w) >= 0) >= gamma];
P = optimizer(Model,wmean'*wmean-gamma.^2,sdpsettings('solver','fmincon'),gamma,wmean'*wmean-gamma^2)


w = sdpvar(2,1);
wmean = [3;4];
a = [1;3];
sdpvar s

Model = [uncertain(w,'normal',wmean,1), probability(a'*w >= s) >= .5,
                                        probability(a'*w <= s) >= .5];
optimize(Model)

sdpvar s1 s2
Model = [uncertain(w,'exponential',wmean), probability(a'*w >= s1) >= .5,
                                           probability(a'*w <= s2) >= .5];
optimize(Model,abs(s1-s2))                                   

sdpvar s1 s2
Model = [uncertain(w,@doublesidedexp,wmean), probability(a'*w >= s1) >= .5,
                                           probability(a'*w <= s2) >= .5];
optimize(Model,abs(s1-s2))                                   
      

                                    
w = sdpvar(2,1);
wmean = [3;4];
a = [1;3];
sdpvar s

Model = [uncertain(w,'normal',wmean,1), probability(a'*w >= s) >= .5,
                                        probability(a'*w <= s) >= .5];                                    

                
% Optimal sensor fusion of normal
clear all
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'normalm',[0;0;0;0],diag([1;2;3;4])),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))


% Optimal sensor fusion of normal, merge distributions automatically
yalmip('clear')
x = sdpvar(4,1);
w = sdpvar(2,1);
q = sdpvar(2,1);
sdpvar s m
Model = [   uncertain(w(1),'normal',0,[1]),        
            uncertain(w(2),'normal',0,[2]),                                             
            uncertain(q ,'normal',[0;0],[3;4]),
            probability(x'*[w;q] >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))

% Optimal sensor fusion of normal, merge distributions automatically
yalmip('clear')
x = sdpvar(4,1);
w = sdpvar(2,1);
q = sdpvar(2,1);
sdpvar s m1 m2 m3 m4
Model = [uncertain(w(1),'normal',m1,1),
         uncertain(w(2),'normal',m2,1),         
         uncertain(q ,'normal',[m3;m4],[1;1]),         
         uncertain(m1,'normal',0,.1),
         uncertain(m2,'normal',0,.1),
         uncertain(m3,'normal',0,.1),
         uncertain(m4,'normal',0,.1),
         probability(x'*[w;q] >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))


% Optimal sensor fusion of normal
clear all
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'normalm',[0;0;0;0],diag([1;2;3;4])),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))

clear all
x = sdpvar(4,1);
w1 = sdpvar(4,1);
w2 = sdpvar(4,1);

sdpvar s
Model = [uncertain(w1,'exponential',([1;2;3;4])),
         uncertain(w2,'exponential',([1;2;3;4])),
         probability(x'*(w1-w2) >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))

sdpvar s
Model = [uncertain(w1,@doublesidedexp,([1;2;3;4])),         
         probability(x'*(w1) >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))

% Optimal sensor fusion of normal, factorised covariance specified
clear all
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'normalf',[0;0;0;0],sqrtm(diag([1;2;3;4]))),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''))

% Optimal sensor fusion of unkown with specified normal
clear all
x = sdpvar(4,1);
w = sdpvar(4,1);
sdpvar s
Model = [uncertain(w,'moment',[0;0;0;0],diag([1;1;1;1])),
         probability(x'*w >= -1) >= s, sum(x)==1, x>=0];
optimize(Model,-s,sdpsettings('debug',1,'solver',''));