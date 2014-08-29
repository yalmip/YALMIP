function [x,w,x_variables,w_variables,aux_variables,F,failure] = robust_classify_variables(F,h,ops,w);

failure = 0;
% Variables before expanding nonlinear operators
initial_variables = unique([depends(F) depends(h)]);

% Uncertain variables
w_variables = getvariables(w);

% Decision variables
x_variables = 1:yalmip('nvars');
x_variables = setdiff(x_variables,w_variables);

% Any fancy modelling goin on here? If so, we need to expand the model, and
% associate the new lifted variables either to the uncertainty set or to
% the set of general variables. There are some classes of variables
% 1. Uncertain variables
% 2. Original decision variables
% 4. Lifted variables used for uncertain variables in uncertainty description
% 5. Lifted variables used for uncertain variables in uncertain model
% 6. Lifted variables used for mixed  variables in uncertain model
% FIX : This code is currently a miserable hack and should be cleaned up
extended = yalmip('extvariables');
aux_variables = [];
if ~isempty(extended) & any(ismember(initial_variables,extended))

    % Remove auxilliary variables (defined by YALMIP)
    aux_variables = intersect(extended,x_variables);
    % These variables are original decision variables
    x_variables = setdiff(x_variables,aux_variables);

    % Original constraint index
    id = getlmiid(F(find(~lifted(F))));
    
    % This is a tweak to allow epxansion of bilinear terms in robust problems,
    % is expression such as abs(x*w) < 1 for all -1 < w < 1
    ops.expandbilinear = 1;

    % Expand the model to model norms etc
    [F,failure,cause] = expandmodel(F,h,ops,w);
    if failure % Convexity propgation failed
        interfacedata = [];
        recoverdata = [];
        solver = '';
        diagnostic.solvertime = 0;
        diagnostic.problem = 14;
        diagnostic.info = yalmiperror(14,cause);
        x = [];
        w = [];
        x_variables=[];
        w_variables=[];
        F = [];
        return
    end

    % Auxillary variables (introduced by YALMIP). These variables are lifting
    % variables that are used to model advanced constructs suck as norms etc in
    % the nonlinear operator framework
    new_aux_variables = setdiff(unique([depends(F) depends(h)]),initial_variables);
    aux_variables = setdiff(union(aux_variables,new_aux_variables),w_variables);

    % We know try to find the auxillary variables that are used to model
    % functions that only depends on w
    F_lifted = F(find(~ismember(getlmiid(F),id)));
    aux_variables_w = aux_variables;
    for i = 1:length(F_lifted)
        variables = depends(F_lifted(i));
        with_x = any(ismember(variables,x_variables));
        with_w = any(ismember(variables,w_variables));
        if with_x
            aux_variables_w = setdiff(aux_variables_w,intersect(variables,aux_variables));
            x_variables = union(x_variables,intersect(variables,aux_variables));
        end
    end

    aux_variables_x = setdiff(aux_variables,aux_variables_w);
    x_variables = union(aux_variables_x, x_variables);
    w_variables = union(aux_variables_w, w_variables);
end

x_variables = intersect([depends(F) depends(h)],x_variables);
w_variables = intersect([depends(F) depends(h)],w_variables);

x = recover(x_variables);
x = recover(setdiff(depends(x),w_variables));
x_variables = getvariables(x);

% % Could be some of these two cases
% % F = (abs(x) + w < 1) + (norm(w,1) < 0.3);
% % F = (x + abs(w) < 1) + (-2 < w < 2) ;
% 
% w_variables = setdiff(w_variables,heh);
% aux_variables = union(aux_variables,heh);

w = recover(w_variables);