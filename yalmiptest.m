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
            clc
            disp('You seem to have some other toolbox with a function called constraint.m');
            disp('Delete that toolbox, or delete the function/class, or change path so that YALMIP is on top.');
            disp(detected{1})
        end
        return
    end
end

detected = which('yalmip.m','-all');
% Will not work in Octave as Octave only reports first item found?
if isa(detected,'cell')
    if length(detected)>1
        clc
        disp('You seem to have multiple installations of YALMIP in your path.')
        disp('Please correct this...');
        disp(detected)
        return
    end
end

% Pagination really doesn't work well with solvers
more off

if exist('OCTAVE_VERSION', 'builtin') 
    OctaveRunning = 1;
else
    OctaveRunning = 0;
end

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
    clc
    disp('The directory yalmip/solvers is not in your path.')
    disp('These must be in path:')
    disp('                      yalmip/');
    disp('                      yalmip/extras');
    disp('                      yalmip/operators');
    disp('                      yalmip/modules');
    disp('                      yalmip/solvers');
    disp('See <a href="https://yalmip.github.io/tutorial/installation/">installation guide</a>')
    return
end

foundstring = {'not found','found','internal'};
teststring = {'-failed','+passed'};
if ~donttest
    header = {'Solver','Version','Status','Unit test'};
else
    header = {'Solver','Version','Status'};
end

[solvers,found] = getavailablesolvers(0);
j = 1;
status = ones(length(solvers),1);
for i = 1:length(solvers)
    if solvers(i).show
        s = {solvers.tag};
        same = find(strcmpi(solvers(i).tag,s));
        % Find all instances of same solver different versions
        for k = setdiff(same,i)
            % No reason to show for all versions
            solvers(k).show = 0;
        end
        
        data{j,1} = upper(solvers(i).tag);
        found_versions = same(find(found(same)));
        if ~isempty(found_versions)
            version = solvers(min(found_versions)).version;
            if ~any(strcmpi(version,{'geometric','standard'}))
                data{j,2} = solvers(min(found_versions)).version;
            end
        end
        status(j) = found(i)+1+solvers(i).builtin;
        data{j,3} = foundstring{found(i)+1+solvers(i).builtin};
        j = j+1;
    end
end

if isa(prefered_solver,'char')
    ops = sdpsettings('Solver',prefered_solver);
else
    ops = prefered_solver;
end
ops.saveyalmipmodel = 1;

if ~((nargin==2) & (ops.verbose==0))
    [sortedName,loc] = sort({data{:,1}});    
    loc = [loc(find(status(loc)==3)) loc(find(status(loc)==2)) loc(find(status(loc)==1))];    
    dataSorted = reshape({data{loc,:}},[],3);    
    yalmiptable({'Searching for installed solvers'},header,dataSorted);
    disp(' ')
end
if nargin<2
    disp('Press any key to continue test')
    pause
end

i=1;
test{i}.fcn  = 'test_core';
test{i}.desc = 'Core functionalities';
i = i+1;

test{i}.fcn  = 'test_linear_programming';
test{i}.desc = 'Linear programming (LP)';
i = i+1;

test{i}.fcn  = 'test_quadratic_programming';
test{i}.desc = 'Quadratic programming (QP)';
i = i+1;

test{i}.fcn  = 'test_socp_programming';
test{i}.desc = 'Second-order cone programming (SOCP)';
i = i+1;

test{i}.fcn  = 'test_semidefinite_programming';
test{i}.desc = 'Semidefinite programming (SDP)';
i = i+1;

test{i}.fcn  = 'test_geometric_programming';
test{i}.desc = 'Geometric programming (GP)';
i = i+1;

test{i}.fcn  = 'test_nonlinear_programming';
test{i}.desc = 'Nonlinear programming (NLP)';
i = i+1;

test{i}.fcn  = 'test_nonlinear_semidefinite_programming';
test{i}.desc = 'Nonlinear SDP (NLSDP)';
i = i+1;

test{i}.fcn  = 'test_exponential_cone_programming';
test{i}.desc = 'Exponential cone programming (ECP)';
i = i+1;

test{i}.fcn  = 'test_milinear_programming';
test{i}.desc = 'Mixed-integer LP (MIQP)';
i = i+1;

test{i}.fcn  = 'test_miquadratic_programming';
test{i}.desc = 'Mixed-integer QP (MIQP)';
i = i+1;

test{i}.fcn  = 'test_misocp_programming';
test{i}.desc = 'Mixed-integer SOCP (MISOCP)';
i = i+1;

test{i}.fcn  = 'test_nonconvex_quadratic_programming';
test{i}.desc = 'Global nonconvex quadratic programming';
i = i+1;

test{i}.fcn  = 'test_nonconvex_global_programming';
test{i}.desc = 'Global nonconvex programming';
i = i+1;

pass_strings = {'Error','Passed','Solver not available'};

% Run test-problems
for i = 1:length(test)
    try
        if ops.verbose
            disp(' ');
            disp(['Testing function ' test{i}.fcn]);
            disp(' ');
        end
        % First make call to figure out solver
        info = eval([test{i}.fcn '(ops)']);
        if ~OctaveRunning
            sols{i} = addLink(upper(cleanversion(info.yalmipmodel.solver.tag)));
        else
            sols{i} = cleanversion(info.yalmipmodel.solver.tag);
        end               
        % And now run all the tests
        info = eval([test{i}.fcn '(ops)']);
        pass(i) = 1;
        results{i}='Success';
    catch
        pass(i) = 0;
        results{i} = 'Failed';
        sols{i} = '';        
    end
end

clear data;
header = {'Test','Status', 'Solver'};
for i = 1:length(pass)
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
    only_fmincon = strcmpi(m.solver.tag,'fmincon-standard');
else
    only_lmilab = 0;
    only_fmincon = 0;
end
if isempty(m)
    disp('You do not have any LMI solver installed')
    disp(' If you intend to solve LMIs you must install a solver.')
else
    if only_lmilab
        disp('You do not have any good LMI solver installed')
        disp(' (only found <a href="https://yalmip.github.io/solver/lmilab/">LMILAB which should be avoided in YALMIP</a>).')
        disp('If you intend to solve LMIs, please install a better solver.')
    elseif only_fmincon
        disp('You do not have any LMI solver installed')
        disp(' (YALMIP will use a nonlinear solver which cannot be expected to work)')
        disp(' If you intend to solve LMIs you must install a solver.')
    end
end
x = binvar(1);[p,aux1,aux2,m] = export(x>=0,x,[],[],[],0);
if isempty(m)
    disp('You do not have any LP/MILP solver installed')
    disp(' If you intend to solve LPs/MILPs, you have to install one.')
else
    only_bnb = strcmpi(m.solver.tag,'bnb');
    if only_bnb
        disp('You do not have any MILP solver installed')
        disp(' (only found internal <a href="https://yalmip.github.io/solver/bnb/">BNB</a>).')
        disp(' If you intend to solve MILP, please install a better solver.')
    end
end
x = binvar(1);[p,aux1,aux2,m] = export(x>=0,x^2,[],[],[],0);
if isempty(m)
    disp('You do not have any QP/SOCP/MIQP/MISOCP solver installed')
    disp(' If you intend to solve QP/SOCP/MIQP/MISOCP you must install solver.')
else
    only_bnb = strcmpi(m.solver.tag,'bnb');
    if only_bnb
        disp('You do not have any MIQP/MISOCP solver installed (only found internal <a href="https://yalmip.github.io/solver/bnb/">BNB</a>)')
        disp(' If you intend to solve MIQP/MISOCP, please install a better solver.')
    end
end
disp('See <a href="https://yalmip.github.io/allsolvers">guide on interfaced solvers</a>')



function info = test_core(ops)
if nargout > 0
    % Fake
    info.yalmipmodel.solver.tag = '';
    return
end
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


function sol = test_semidefinite_programming(ops)
%%
t = sdpvar(1,1);
Y = sdpvar(2,2);
F = [Y<=t*eye(2), Y>=[1 0.2;0.2 1]];
sol = optimize(F,t,ops);
if nargout > 0
    return
end
assert(ismember(sol.problem,[0 3 4 5]) && abs(value(t)-1.2)<=1e-4);
%%
A = [1 0;0.4 1];
B = [0.4;0.08];
L = [1.9034 1.1501];
Y = sdpvar(2,2);
F = [Y Y*(A-B*L)';(A-B*L)*Y Y]>=0;
F = F+[L*Y*L'<=1];
sol = optimize(F,-logdet(Y),ops);
Y = double(Y);
assert(ismember(sol.problem,[0 3 4 5]) && norm(value(Y)-[2.9957 -4.1514;-4.1514 6.2918])<=1e-3)
%
x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);
X = [[x 1 2];[1 y 3];[2 3 100]];
F = [X>=0,x>=10,y>=0,z>=0, x<=1000, y<=1000,z<=1000];
sol = optimize(F,x+y+z,ops);
x   = double(x);
y   = double(y);
z   = double(z);
assert(ismember(sol.problem,[0 3 4 5]) && norm(value([x;y;z])-[10;0.1787;0])<=1e-3)
%%
x = sdpvar(1,1);
z = sdpvar(1,1);
X = [[x 2];[2 z]];
F = [X>=0, x>=0,z>=0,x<=10,z<=10];
sol = optimize(F,x+z,ops);
assert(~sol.problem)
%%
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
assert(sol.problem == 0 && abs(value(t)-14.8361)<=1e-3)
%%
n = 5;
P = magic(n);
Z = sdpvar(n,n,'toeplitz');
t = sdpvar(n,n,'full');
F = (P(:)-Z(:)<=t(:))+(P(:)-Z(:)>=-t(:));
sol = optimize(F,sum(sum(t)),ops);
assert(~sol.problem)
assert(abs(value(sum(sum(t)))-156)<=1e-3)
%%
n = 5;
P = magic(n);
Z = sdpvar(n,n,'toeplitz');
t = sdpvar(n,n,'full');
resid = P-Z;resid = resid(:);
sol = optimize([],resid'*resid,ops);
assert(~sol.problem)
assert(abs(value(resid'*resid)-1300)<=1e-2)
%%
A = magic(6);
A = A*A';
P = sdpvar(6,6);
sol = optimize((A'*P+P*A <= -P) + (P>=eye(6)),trace(P),ops);
assert(sol.problem==1)

function sol = test_linear_programming(ops)
%%
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
if nargout > 0
    return
end
assert(~sol.problem)
assert(abs(value(sum(t))-12.66666)<=1e-3)
%%
t = sdpvar(1,1);
F = [t>=0, t<=-10];
sol = optimize(F,t,ops);
assert(sol.problem~=0)

function sol = test_socp_programming(ops)
x = sdpvar(2,1);
a = [0;1];
b = [1;1];
F = norm(x-a)<=1;
F = F+[norm(x-b) <= 1];
sol = optimize(F,sum(x),ops);
if nargout > 0
    return
end
assert(~sol.problem)
assert(abs(value(sum(x))-0.58578)<=1e-3)
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
assert(~sol.problem)
x = double(x);
y = double(y);
z = double(z);
assert(abs(sum(x)-0.27592)<=1e-3);
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
assert(~sol.problem)
x = double(x);
assert(abs(sum(x)-0.59)<=1e-3)

function sol = test_misdp_programming(ops)
x = intvar(4,1);
e = magic(4)*x-1;
sdpvar t
obj = t;
sol = optimize([t e';e eye(4)]>=0,obj,ops);
if nargout > 0
    return
end
assert(~sol.problem)
assert(abs(value(obj)-4)<=1e-3)

function sol = test_misocp_programming(ops)
x = intvar(4,1);
obj = norm(magic(4)*x-1,2);
sol = optimize([-5 <= x <= 5],obj,ops);
if nargout > 0
    return
end
assert(~sol.problem)
assert(abs(value(obj)-2)<=1e-3)

function sol = test_miquadratic_programming(ops)
x = intvar(4,1);
obj = norm(magic(4)*x-1,2)^2
sol = optimize([-5 <= x <= 5],obj,ops);
if nargout > 0
    return
end
assert(~sol.problem)
assert(abs(value(obj)-4)<=1e-3)

function sol = test_milinear_programming(ops)
x = intvar(4,1);
obj = norm(magic(4)*x-1,1)
sol = optimize([-5 <= x <= 5],obj,ops);
if nargout > 0
    return
end
assert(~sol.problem)
assert(abs(value(obj)-2)<=1e-3)

function sol = test_quadratic_programming(ops)
x = sdpvar(10,1);
sol = optimize([sum(x)==2, -1 <= x <= 1],x'*x,ops);
if nargout > 0
    return
end
assert(~sol.problem)
assert(abs(value(obj)-.4)<=1e-3)

function sol = test_nonconvex_quadratic_programming(ops)
x = sdpvar(10,1);
ops.forceglobal = 1;
sol = optimize([sum(x)==2, -1 <= x <= 1],-x'*x,ops);
if nargout > 0
    return
end
assert(~sol.problem)
assert(abs(value(obj)-10)<=1e-3)

function sol = test_nonconvex_global_programming(ops)
x = sdpvar(3,1);
ops.forceglobal = 1;
sol = optimize([sum(x.^3)==2, -1 <= x <= 1],-x'*x,ops);
if nargout > 0
    return
end
assert(~sol.problem)
assert(abs(value(obj)-2.28)<=1e-1)


function sol = test_nonlinear_semidefinite_programming(ops)
%%
A = [-1 2;-3 -4];
P = sdpvar(2,2);
alpha = sdpvar(1,1);
F = (P>=eye(2))+(A'*P+P*A <= -2*alpha*P)+(alpha >= 0);
sol = optimize([F,P(:) <= 100],-alpha,ops);
if nargout > 0
    return
end
assert(~sol.problem);
assert(abs(value(alpha)-2.5)<= 1e-3)

function sol = test_geometric_programming(ops)
%%
t1 = sdpvar(1,1);
t2 = sdpvar(1,1);
t3 = sdpvar(1,1);
t = [t1 t2 t3];
obj = (40*t1^-1*t2^-0.5*t3^-1)+(20*t1*t3)+(40*t1*t2*t3);
F = ((1/3)*t1^-2*t2^-2+(4/3)*t2^0.5*t3^-1 <= 1);
F = [F, t>=0];
sol = optimize(F,obj,ops);
if nargout > 0
    return
end
assert(~sol.problem)
assert(abs(value(D)-133.1371)<=1e-3)

function sol = test_nonlinear_programming(ops)
sdpvar x y
sol = optimize(x^2 + x^4 + exp(x) <= 1, x^2+y^2,ops)
if nargout > 0
    return
end
assert(sol.problem == 0);

function sol = test_exponential_cone_programming(ops)
sdpvar x y z
sol = optimize([expcone([x;2;z]),x==1],z,ops)
if nargout > 0
    return
end
assert(sol.problem == 0);

function html = addLink(x)

if length(x)>0
    html = ['<a href="https://yalmip.github.io/solver/' lower(x) '">' upper(x) '</a>'];  
else
    html = '';
end

function x = cleanversion(x)

s = strfind(x,'-');
if ~isempty(s)
    x=x(1:s-1);
end