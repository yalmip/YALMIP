function output = callpenlab(model)

model = yalmip2nonlinearsolver(model);

if ~model.derivative_available
    disp('Derivate-free call to penlab not supported')
    error('Derivate-free call to penlab not implemented')
end
if model.options.savedebug
    save penlabdebug model
end

penm = [];
penm.Nx=length(model.linearindicies);


if length(nnz(model.K.s)>0)
    top = 1 + model.K.f + model.K.l + sum(model.K.q);
    for i = 1:length(model.K.s)
        model.vecF{i} = model.F_struc(top:top+model.K.s(i)^2-1,:);
        top = top + model.K.s(i)^2;
    end
end

showprogress('Calling PENLAB',model.options.showprogress);

% These are needed to avoid recomputation due to penlab double call to get
% f and df, and g and dg
global latest_x_f
global latest_x_g
global latest_df
global latest_f
global latest_G
global latest_g
global latest_xevaled
global latest_x_xevaled
latest_G = [];
latest_g = [];
latest_x_f = [];
latest_x_g = [];
latest_xevaled = [];
latest_x_xevaled = [];

penm.userdata=model;

penm.objfun =  @(x,Y,userdata) penlab_callback_f(x,userdata);
penm.objgrad = @(x,Y,userdata) penlab_callback_df(x,userdata);
penm.objhess = @(x,Y,userdata) penlab_callback_df2(x,userdata);

if length(model.b)>0 | length(model.beq)>0
    penm.lbg = [model.beq;model.b-inf];
    penm.ubg = [model.beq;model.b];
    penm.NgLIN = length(model.beq) + length(model.b);
    penm.confun =  @(x,Y,userdata) penlab_callback_con(x,userdata);
    penm.congrad = @(x,Y,userdata) penlab_callback_dcon(x,userdata);
    penm.conhess = @(x,Y,k,userdata) penlab_callback_dcon2(x,k,userdata);
end

if model.K.s(1)>0
    penm.NALIN=length(model.K.s);
    penm.lbA=zeros(find(model.K.s)>0,1);
    penm.mconfun  = @(x,Y,k,userdata) penlab_callback_matrixG(x, k,userdata);
    penm.mcongrad = @(x,Y,k,i,userdata) penlab_callback_matrixdG(x,k,i,userdata);
end

prob=penlab(penm);  

solvertime = tic;
prob.solve();       
solvertime = toc(solvertime);
x = prob.x;

% Duals currently not supported
lambda = [];

switch 0
    case {0,1}
        problem = 0;
    case {2}
        problem = 1;
    case {-1}
        problem = 3;
    case {3,4,-2,-3}
        problem = 4;
    case {-11,-12,-13}
        problem = 7;
    case {-10,-100,-101,-102,-199}
        problem = 11;
    otherwise
        problem = -1;
end

% Internal format for duals
D_struc = [];

% Save all data sent to solver?
if model.options.savesolverinput
    solverinput.model = model;
else
    solverinput = [];
end

% Save all data from the solver?
if model.options.savesolveroutput
    solveroutput.x = xout;  
    solveroutput.info = info;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'PENLAB',solverinput,solveroutput,solvertime);


