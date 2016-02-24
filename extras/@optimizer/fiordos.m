function solver = fiordos(self)

% The first dimin equalities are used to fix some parameters
aux = self.model.F_struc(1:self.dimin,:);

% remove these 
self.model.F_struc(1:self.dimin(1),:) = [];
self.model.K.f = self.model.K.f - self.dimin(1);

% Extract all bounds and move from Ab to ub/lb
temp =  presolve_bounds_from_modelbounds(self.model,1);   
self.model = temp;
% Any bound constraints?
n = length(temp.lb);
bounded_ub = find(~isinf(temp.ub));
bounded_lb = find(~isinf(temp.lb));
if length(bounded_lb) == n && length(bounded_ub) == n
    % All have finite bounds.
    X = EssBox(n, 'l',temp.lb, 'u',temp.ub);
    X = SimpleSet(X);   
elseif  ~isempty(bounded_ub) || ~isempty(bounded_lb)
    % Only some bounds. Keep them as general constraints
    X = [];
    for i = 1:n
        if isinf(temp.lb(i)) && isinf(temp.ub(i))
            X{end+1} = EssRn(1);
        elseif  ~isinf(temp.lb(i)) && ~isinf(temp.ub(i))
            X{end+1} =  EssBox(1, 'l',temp.lb(i), 'u',temp.ub(i));
        elseif isinf(temp.ub(i))            
            X{end+1} =  EssRnplus(1, 'shift', temp.lb(i));
        else
            X{end+1} =  EssRnplus(1, 'shift', temp.ub(i),'rot',-1);
        end
    end
    X = SimpleSet(X{:});
else
    % No bounds
    X = EssRn(n);
    X = SimpleSet(X);
end

% Put back a placeholder for fixing the parameters
self.model.F_struc = [aux;self.model.F_struc];
self.model.F_struc(1:self.dimin(1),1)=0; % Now fixed to zero
self.model.K.f = self.model.K.f + self.dimin(1);

op = OptProb('H',2*full(self.model.Q), 'g',self.model.c, 'X',X, 'Ae',-self.model.F_struc(1:self.model.K.f,2:end), 'be','param'); 

%instantiate solver
s = Solver(op,'approach','primal-dual');
%, 'algoOuter','gm', 'algoInner','fgm'); 
%optionally change settings, e.g.
%-maximum number of iterations
%s.setSettings('algoOuter', 'maxit',10000);
%s.setSettings('algoInner', 'maxit',9000);
%-gradient-map stopping criterion
%s.setSettings('algoOuter', 'stopg',true, 'stopgEps',1e-3); 
%s.setSettings('algoInner', 'stopg',true, 'stopgEps',1e-5);
%generate solver code  
s.generateCode('prefix','demo_','forceOverwrite',1);
demo_mex_make();

compiledsolver = @(x)demo_mex(x);
b0 = self.model.F_struc(1:self.model.K.f,1);
B0 = [eye(self.dimin(1));zeros(self.model.K.f-self.dimin(1),self.dimin(1))];
map = self.map;
dimout = self.dimout;
mask = self.mask{1};
solver = @(x)(fiordos_call(compiledsolver,x,B0,b0,mask,map,dimout));
