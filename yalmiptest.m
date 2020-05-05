function out = yalmiptest(prefered_solver,auto)
%YALMIPTEST Runs a number of test problems.
%
%   YALMIPTEST is recommended when a new solver or a new version
%   of YALMIP installed.
%
%   EXAMPLES
%    YALMIPTEST               % Without argument, default solver used
%    YALMIPTEST('solver tag') % Test with specified solver
%    YALMIPTEST(options)      % Test with specific options structure from
%
%   See also SDPSETTINGS

if ~exist('sedumi2pen.m')
    disp('Add /yalmip/extras etc to your path first...')
    disp('Read the <a href="https://yalmip.github.io/tutorial/installation/">Installation notes</a>.')    
    return
end

if ~exist('callsedumi.m')
    disp('Still missing paths...Just do an addpath(genpath(''yalmiprootdirectory''));')
    disp('Read the <a href="https://yalmip.github.io/tutorial/installation/">Installation notes</a>.')        
    return
end

% SDPT3 has a function called constraint.m which causes issues
detected = which('constraint.m');
if isa(detected,'cell') 
    if length(detected)>0              
        if isempty(strfind(detected{1},'extras\@constraint'))
            disp('You seem to have some other toolbox with a function called constraint.m');
            disp('Delete that toolbox, or delete the function/class, or change path so that YALMIP is on top.');
            detected{1}            
        end
        return
    end
end

detected = which('yalmip.m','-all');
% Will not work in Octave as Octave only reports first item found?
if isa(detected,'cell') 
    if length(detected)>1
        disp('You seem to have multiple installations of YALMIP in your path. Please correct this...');
        detected
        return
    end
end

% Pagination really doesn't work well with solvers 
more off

donttest = 0;
if (nargin==1) && isa(prefered_solver,'char') && strcmp(prefered_solver,'test')
    donttest = 0;
    prefered_solver = '';    
else
    donttest = 1;
end

if nargin==0
    prefered_solver = '';
else
    if ~(isa(prefered_solver,'struct') | isa(prefered_solver,'char'))
        error('Argument should be a solver tag, or a sdpsettings structure');
    end
    if isa(prefered_solver,'char') 
        donttest = 1;
    end
end

if ~(exist('callsedumi')==2)
    disp('The directory yalmip/solvers is not in your path.')
    disp('Put yalmip/, yalmip/solvers, yalmip/extras and yalmip/demos in your MATLAB path.');
    return
end

foundstring = {'not found','found'};
teststring = {'-failed','+passed'};
if ~donttest
    header = {'Solver','Version/module','Status','Unit test'};
else
    header = {'Solver','Version/module','Status'};
end    

[solvers,found] = getavailablesolvers(0);
solvers = solvers([find(found);find(~found)]);
found = [found(find(found));found(find(~found))];
j = 1;
for i = 1:length(solvers)
    if solvers(i).show
        data{j,1} = upper(solvers(i).tag); 
        data{j,2} = solvers(i).version;
        if length(solvers(i).subversion)>0
            data{j,2} = [data{j,2} ' ' solvers(i).subversion];
        end
        data{j,3} = foundstring{found(i)+1};   
        if ~donttest
            if found(i)
                if options.verbose
                    disp(['Testing ' solvers(i).tag '...']);    
                end
                try
                    if solvers(i).maxdet
                        pass = lyapell(sdpsettings('solver',solvers(i).tag,'verbose',0));            
                    else
                        if solvers(i).sdp
                            pass = stabtest(sdpsettings('solver',solvers(i).tag,'verbose',0));
                        else
                            pass = feasiblelp(sdpsettings('solver',solvers(i).tag,'verbose',0));
                        end
                    end
                    data{j,4} = teststring{pass+1};
                catch
                    data{j,4} = '-failed';
                end
            else
                data{j,4} = 'not tested';
            end
        end
        j = j+1;
    end
end

if isa(prefered_solver,'char')
    ops = sdpsettings('Solver',prefered_solver);
else
    ops = prefered_solver;
end

if ~((nargin==2) & (ops.verbose==0))
    [sortedName,loc] = sort({data{:,1}});
    dataSorted = reshape({data{loc,:}},[],3);
    yalmiptable({'Searching for installed solvers'},header,dataSorted);
    disp(' ')
end
if nargin<2
    disp('Press any key to continue test')
    pause
end

i=1;
test{i}.fcn  = 'testsdpvar';
test{i}.desc = 'Core functionalities';
i = i+1;

test{i}.fcn  = 'feasiblelp'; 
test{i}.desc = 'LP';
i = i+1;

test{i}.fcn  = 'toepapprox'; 
test{i}.desc = 'LP';
i = i+1;

test{i}.fcn  = 'feasibleqp'; 
test{i}.desc = 'QP';
i = i+1;

test{i}.fcn  = 'toepapprox2'; 
test{i}.desc = 'QP';
i = i+1;


test{i}.fcn  = 'socptest1'; 
test{i}.desc = 'SOCP';
i = i+1;

test{i}.fcn  = 'socptest2'; 
test{i}.desc = 'SOCP'; 
i = i+1;

test{i}.fcn  = 'socptest3'; 
test{i}.desc = 'SOCP';   
i = i+1;

test{i}.fcn  = 'complete';
test{i}.desc = 'SDP';
i = i+1;

test{i}.fcn  = 'complete_2'; 
test{i}.desc = 'SDP';
i = i+1;

test{i}.fcn  = 'maxcut';
test{i}.desc = 'SDP';
i = i+1;

test{i}.fcn  = 'feasible'; 
test{i}.desc = 'SDP';
i = i+1;

test{i}.fcn  = 'lyapell'; 
test{i}.desc = 'MAXDET';
i = i+1;

test{i}.fcn  = 'lyapell2'; 
test{i}.desc = 'MAXDET';
i = i+1;

%test{i}.fcn  = 'circuit1'; 
%test{i}.desc = 'GP';
%i = i+1;

test{i}.fcn  = 'infeasible'; 
test{i}.desc = 'Infeasible LP';
i = i+1;

test{i}.fcn  = 'infeasibleqp'; 
test{i}.desc = 'Infeasible QP';
i = i+1;

test{i}.fcn  = 'infeasiblesdp'; 
test{i}.desc = 'Infeasible SDP';
i = i+1;

test{i}.fcn  = 'momenttest'; 
test{i}.desc = 'Moment relaxation';
i = i+1;

test{i}.fcn  = 'sostest'; 
test{i}.desc = 'Sum-of-squares';
i = i+1;

test{i}.fcn  = 'bmitest'; 
test{i}.desc = 'Bilinear SDP';
i = i+1;



pass_strings = {'Error','Passed','Solver not available'};

tt = cputime;

% Run test-problems
for i = 1:length(test)
    try
        t=cputime;           
        if ops.verbose
            disp(' ');
            disp(['Testing function ' test{i}.fcn]);
            disp(' ');
        end       
        [pp,ss,res] = eval([test{i}.fcn '(ops)']);
        pass(i) = pp;
        sols{i} = ss.info;
        results{i}=res;
        ttime(i) = cputime-t;
    catch
        pass(i) = 0;   
        results{i} = 'NAN';
        sols{i} = 'Unknown problem in YALMIP';
        ttime(i) = cputime-tt;
    end
end
totaltime = cputime-tt;

clear data;
header = {'Test','Solution', 'Solver message'};
for i = 1:length(pass)
    thetime =  num2str(ttime(i),4);
    data{i,1} = test{i}.desc;
    data{i,2} = results{i};
    data{i,3} = sols{i};
end
if ops.verbose
    disp(' ');
end
formats{1}.data.just = 'right';
formats{2}.data.just = 'right';
formats{3}.data.just = 'right';

formats{1}.header.just = 'right';
formats{2}.header.just = 'right';
formats{3}.header.just = 'right';

clc
yalmiptable([],header,data,formats)

% Test if any LMI solver is installed.
x = sdpvar(2);[p,aux1,aux2,m] = export(x>=0,[],[],[],[],0);
if ~isempty(m)
  only_lmilab = strcmpi(m.solver.tag,'lmilab');
else
  only_lmilab = 0;
end
x = binvar(1);[p,aux1,aux2,m] = export(x>=0,[],[],[],[],0);
if ~isempty(m)
  only_bnb = strcmpi(m.solver.tag,'bnb');
else
  only_bnb = 0;
end
if only_lmilab 
 disp('You do not have any efficient LMI solver installed (only found <a href=" http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Solvers.LMILAB">LMILAB</a>).')
 disp('If you intend to solve LMIs, please install a better solver.')
end
if only_bnb 
 disp('You do not have any efficient MILP solver installed (only found internal <a href=" http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Solvers.BNB">BNB</a>).')
 disp('If you intend to solve MILPs, please install a better solver.')
end
if only_lmilab  || only_bnb
  disp('See <a href=" http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Solvers">Interfaced solvers in YALMIP</a>')
end


function [pass,sol,result] = testsdpvar(ops)

% Test the sdpvar implementation
pass = 1;
sol.info = yalmiperror(0,'YALMIP');
try
    x = sdpvar(2,2);
    x = sdpvar(2,2,'symmetric');
    x = sdpvar(2,2,'full');
    x = sdpvar(2,2,'toeplitz');
    x = sdpvar(2,2,'hankel');
    x = sdpvar(2,2,'skew');
    if ~ishermitian(sdpvar(2,2,'hermitian','complex'))
        error('bug')
    end
    if ~issymmetric(sdpvar(2,2,'symmetric','complex'))
        error('bug')
    end
    if ~isreal(real(sdpvar(2,2,'symmetric','complex')))
        error('bug')
    end
    if isreal(sqrt(-1)*real(sdpvar(2,2,'symmetric','complex')))
        error('bug')
    end
    x = sdpvar(2,1,'','co');
    if ~isreal(x'*x)
        error('bug')
    end
    x = sdpvar(2,2,'','co');
    if ~isreal(diag(x'*x))
        error('bug')
    end
    x = sdpvar(1,1);
    y = sdpvar(2,2);
    x*eye(2);
    eye(2)*x;
    y*3;
    3*y;
    x = sdpvar(2,3);
    y = sdpvar(2,3);
    assign(x,randn(2,3));
    z = replace(x,x(1,1:2),[8 9]);
    z = x+y;
    z = x-y;
    z = x+1;
    z = x-1;
    z = x+ones(2,3);
    z = x-ones(2,3);
    z = ones(2,3)-x;
    z = ones(2,3)-x;
    z = eye(2)*x;
    z = x*eye(3);
    z = diag(x);
    z = trace(x(1:2,1:2));
    z = diff(x);
    z = fliplr(x);
    z = flipud(x);
    z = kron(x,eye(3));
    z = kron(eye(3),x);
    z = rot90(x);
    z = sum(x);
    z = diff(x);
    z = x';
    z = x.';
    z = tril(x);
    z = triu(x);
    z = [x y];
    z = [x;y];
    sdpvar x y
    diag([x y])*[x^-1;y^-1];
    assert(isequal([x x;x x]*x-[x x;x x].*x,zeros(2)))
    assert(isequal(trace([x x;x x]*[x y;y x])-(x*x+x*y+y*x+x*x),0))
    
    % Regression ??
    yalmip('clear')
    sdpvar x

    (1+x+x^4)*(1-x^2); 
    
    % Regression complex multiplcation
    A = randn(10,5)+sqrt(-1)*randn(10,5);
    b = randn(10,1)+sqrt(-1)*randn(10,1);
    x = sdpvar(5,1);
    res = A*x-b;
    assert(nnz(clean([res res]'*[res res]-res'*res,1e-8))==0)
    assert(isreal(clean(res'*res,1e-8)))

    assert(isreal(x*x'))
              
    result = 'N/A';
catch
    sol.info = 'Problems';
    result = 'N/A';
    pass = 0;
end


function [pass,sol,result] = feasible(ops) 
t = sdpvar(1,1);
Y = sdpvar(2,2);
F = [Y<=t*eye(2), Y>=[1 0.2;0.2 1]];
sol = optimize(F,t,ops);
pass = ismember(sol.problem,[0 3 4 5]);
if pass
    result = resultstring(t,1.2);
else
    result = 'N/A';
end
    


function [pass,sol,result] = infeasible(ops)
t = sdpvar(1,1);
F = [t>=0, t<=-10];
sol = optimize(F,t,ops);
pass = ~(sol.problem==0);
result = 'N/A';

function [pass,sol,result] = lyapell(ops)
A = [1 0;0.4 1];
B = [0.4;0.08]; 
L = [1.9034 1.1501];

Y = sdpvar(2,2);
F = [Y Y*(A-B*L)';(A-B*L)*Y Y]>=0;
F = F+[L*Y*L'<=1];
sol = optimize(F,-logdet(Y),ops);
Y = double(Y);
pass = ismember(sol.problem,[0 3 4 5]);
if pass
    result = resultstring(Y,[2.9957 -4.1514;-4.1514 6.2918]);
else
    result = 'N/A';
end
%pass = pass & (sum(sum(abs(Y-[2.9957 -4.15;-4.15 6.29])))<0.01);

function [pass,sol,result] = lyapell2(ops)
A = [1 0;0.4 1];
B = [0.4;0.08]; 
L = [1.9034 1.1501];  
Y = sdpvar(2,2);
F = [Y Y*(A-B*L)';(A-B*L)*Y Y]>=0;
F = F+[L*Y*L'<=1];
sol = optimize(F,-logdet(Y),ops);
Y = double(Y);
pass = ismember(sol.problem,[0 3 4 5]);
if pass
    result = resultstring(Y,[2.9957 -4.1514;-4.1514 6.2918]);
else
    result = 'N/A';
end

function [pass,sol,result] = complete(ops)
x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);

X = [[x 1 2];[1 y 3];[2 3 100]];

F = [X>=0,x>=10,y>=0,z>=0, x<=1000, y<=1000,z<=1000];
sol = optimize(F,x+y+z,ops);
x   = double(x);
y   = double(y);
z   = double(z);

pass = ismember(sol.problem,[0 3 4 5]);
result = 'N/A';
if pass
    result = resultstring([x;y;z],[10;0.1787;0]);
else
    result = 'N/A';
end


function [pass,sol,result] = complete_2(ops)
yalmip('clear')
x = sdpvar(1,1);
z = sdpvar(1,1);

X = [[x 2];[2 z]];

F = [X>=0, x>=0,z>=0,x<=10,z<=10];
sol = optimize(F,x+z,ops);
x   = double(x);
z   = double(z);

pass = ismember(sol.problem,[0 3 4 5]);
result = 'N/A';
if pass
    result = resultstring([x;z],[2;2]);
else
    result = 'N/A';
end



function [pass,sol,result]  = maxcut(ops)
% Upper bound on maxcut of a n-cycle
n = 15;
Q = zeros(n);
for i = 1:n-1
    Q(i,i+1) = 1;Q(i+1,i)  = 1;
end
Q(n,1) = 1;Q(1,n) = 1;  
Q = 0.25*(diag(Q*ones(n,1))-Q);

t = sdpvar(1,1);
tau = sdpvar(n,1);

F = t>=0;

M = [[-Q zeros(n,1)];[zeros(1,n) t]];

for i = 1:n
    ei = zeros(n,1);ei(i,1) = 1;
    M = M+tau(i)*[ei*ei' zeros(n,1);zeros(1,n) -1];
end

F = F+[M>=0];
sol = optimize(F,t,ops);

t   = double(t);
tau = double(t);

pass = ismember(sol.problem,[0 3 4 5]);
if pass
    result = resultstring(t,14.8361);
else
    result = 'N/A';
end


function [pass,sol,result] = socptest1(ops)
yalmip('clear')
x = sdpvar(2,1);
a = [0;1];
b = [1;1];
F = norm(x-a)<=1;
F = F+[norm(x-b) <= 1];
sol = optimize(F,sum(x),ops);
pass = ismember(sol.problem,[0 3 4 5]); 

x = double(x);

if pass
    result = resultstring(sum(x),0.58578);
else
    result = 'N/A';
end



function [pass,sol,result] = socptest2(ops)
z = sdpvar(3,1);
x = sdpvar(3,1);
y = sdpvar(3,1);
a = [0;1;0];
b = [1;1;0];
F = norm(x-a)<=1;
F = F+[norm(x-b)<=1];
F = F+[x(1)==0.35];
F = F+[z(2:3)==[5;6]];
sol = optimize(F,sum(x),ops);
pass = ismember(sol.problem,[0 3 4 5]); 

x = double(x);
y = double(y);
z = double(z);

if pass
    result = resultstring(sum(x),0.27592);
else
    result = 'N/A';
end



function [pass,sol,result] = socptest3(ops)
z = sdpvar(2,1);
x = sdpvar(2,1);
y = sdpvar(3,1);
a = [0;1];
b = [1;1];
F = norm(x-a)<=1;
F = F+[norm(x-b)<=1];
F = F+[x(1)==0.35];
F = F+[z(1,end)>=5];
F = F+[z(2,end)<=100];
F = F+[z(2)==5];

sol = optimize(F,sum(x),ops);
pass = ismember(sol.problem,[0 3 4 5]); 

x = double(x);
y = double(y);
z = double(z);
if pass
    result = resultstring(sum(x),0.59);
else
    result = 'N/A';
end



function [pass,sol,result] = feasiblelp(ops)
N = 5;
A = [2 -1;1 0];
B = [1;0];
C = [0.5 0.5];
[H,S] = create_CHS(A,B,C,N);
x = [2;0];
t = sdpvar(2*N,1);
U = sdpvar(N,1);   
Y = H*x+S*U; 
F = (U<=1)+(U>=-1);
F = F+(Y(N)>=-1);  
F = F+(Y(N)<=1); 
F = F+([Y;U]<=t)+([Y;U]>=-t);
sol = optimize(F,sum(t),ops);
pass = ismember(sol.problem,[0 3 4 5]); 
if pass
    result = resultstring(sum(t),12.66666);
else
    result = 'N/A';
end

function [pass,sol,result] = feasibleqp(ops)
N = 5;
A = [2 -1;1 0];
B = [1;0];
C = [0.5 0.5];
[H,S] = create_CHS(A,B,C,N);
x = [2;0];
U = sdpvar(N,1);   
Y = H*x+S*U; 
F = (U<=1)+(U>=-1);
F = F+(Y(N)>=-1);  
F = F+(Y(N)<=1); 
sol = optimize(F,Y'*Y+U'*U,ops);
pass = ismember(sol.problem,[0 3 4 5]); 
if pass
    result = resultstring(Y'*Y+U'*U,26.35248);
else
    result = 'N/A';
end


function [pass,sol,result] = infeasibleqp(ops)
N = 5;
A = [2 -1;1 0];
B = [1;0];
C = [0.5 0.5];
[H,S] = create_CHS(A,B,C,N);
x = [2;0];
U = sdpvar(N,1);   
Y = H*x+S*U; 
F = (U<=1)+(U>=-1);
F = F+(Y(N)>=-1);  
F = F+(Y(N)<=1); 
F = F + (U>=0);
sol = optimize(F,Y'*Y+U'*U,ops);
pass = ismember(sol.problem,[1]); 
result = 'N/A';



function [pass,sol,result] = infeasiblesdp(ops)
A = magic(6);
A = A*A';
P = sdpvar(6,6);
sol = optimize((A'*P+P*A <= -P) + (P>=eye(6)),trace(P),ops); 
pass = (sol.problem==1);
result = 'N/A';

function [pass,sol,result]=toepapprox(ops)

n = 5;
P = magic(n);
Z = sdpvar(n,n,'toeplitz');
t = sdpvar(n,n,'full');
F = (P-Z<=t)+(P-Z>=-t);
sol = optimize(F,sum(sum(t)),ops);
pass = ismember(sol.problem,[0 3 4 5]); 
result = 'N/A';
if pass
    result = resultstring(sum(sum(t)),156);
else
    result = 'N/A';
end


function [pass,sol,result]=toepapprox2(ops)

n = 5;
P = magic(n);
Z = sdpvar(n,n,'toeplitz');
t = sdpvar(n,n,'full');
resid = P-Z;resid = resid(:);
sol = optimize([],resid'*resid,ops);
pass = ismember(sol.problem,[0 3 4 5]); 
result = 'N/A';
if pass
    result = resultstring(resid'*resid,1300);
else
    result = 'N/A';
end

function [pass,sol,result]=momenttest(ops)

x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);

objective = -2*x1+x2-x3;

F = (x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24>=0);
F = F + (4-(x1+x2+x3)>=0);
F = F + (6-(3*x2+x3)>=0);
F = F + (x1>=0);
F = F + (2-x1>=0);
F = F + (x2>=0);
F = F + (x3>=0);
F = F + (3-x3>=0);
sol = solvemoment(F,objective,ops);
pass = ismember(sol.problem,[0 3 4 5]); 
result = 'N/A';
if pass
    result = resultstring(objective,-6);
else
    result = 'N/A';
end

function [pass,sol,result]=sostest(ops)

yalmip('clear')
x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);
F = (sos(1+x^7+x^8+y^4-t));
sol = solvesos(F,-t,ops);
pass = ismember(sol.problem,[0 3 4 5]); 
result = 'N/A';
if pass
    result = resultstring(t,0.9509);
else
    result = 'N/A';
end

function [pass,sol,result]=bmitest(ops)

A = [-1 2;-3 -4];
P = sdpvar(2,2);
alpha = sdpvar(1,1);
F = (P>=eye(2))+(A'*P+P*A <= -2*alpha*P)+(alpha >= 0);
sol = optimize([F,P(:) <= 100],-alpha,ops);
pass = ismember(sol.problem,[0 3 4 5]); 
result = 'N/A';
if pass
    result = resultstring(alpha,2.5);
else
    result = 'N/A';
end


function [pass,sol,result]=circuit1(ops)
x = sdpvar(7,1);

% Data
a     = ones(7,1);
alpha = ones(7,1);
beta  = ones(7,1);
gamma = ones(7,1);
f = [1 0.8 1 0.7 0.7 0.5 0.5]';
e = [1 2 1 1.5 1.5 1 2]';
Cout6 = 10;
Cout7 = 10;

% Model
C = alpha+beta.*x;
A = sum(a.*x);
P = sum(f.*e.*x);
R = gamma./x;

D1 = R(1)*(C(4));
D2 = R(2)*(C(4)+C(5));
D3 = R(3)*(C(5)+C(7));
D4 = R(4)*(C(6)+C(7));
D5 = R(5)*(C(7));
D6 = R(6)*Cout6;
D7 = R(7)*Cout7;

% Constraints
F = (x >= 1) + (P <= 20) + (A <= 100);

% Objective
D = max((D1+D4+D6),(D1+D4+D7),(D2+D4+D6),(D2+D4+D7),(D2+D5+D7),(D3+D5+D6),(D3+D7));

sol = optimize(F,D,ops);

pass = ismember(sol.problem,[0 3 4 5]); 
result = 'N/A';
if pass
    result = resultstring(D,7.8936);
else
    result = 'N/A';
end





function result = resultstring(x,xopt)
if norm(double(x(:))-xopt(:))<=1e-3*(1+norm(xopt(:)))
    result = 'Correct';
else
    result = 'Incorrect';
end

function assert(a)
if ~a
    error('Assertion failed!');
end

