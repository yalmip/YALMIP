function [VariableType,F] = robust_classify_variables_newest(F,h,ops,w);

% Variables before expanding nonlinear operators
initial_variables = unique([depends(F) depends(h)]);

% Uncertain variables
w_variables = getvariables(w);

% Decision variables
x_variables = 1:yalmip('nvars');
x_variables = setdiff(x_variables,w_variables);

% These can have complicated dependence
extended = yalmip('extvariables');

% so we know for sure this are simple variables without any w-dependence
x_variables = setdiff(x_variables,extended);
aux_variables = [];

Graph = [];
aux_without_w_dependence = [];
aux_with_only_w_dependence = [];
aux_with_w_dependence = [];
Graph = [];

if ~isempty(extended) & any(ismember(initial_variables,extended))

    % Remove auxilliary variables (defined by YALMIP)
    aux_variables = intersect(extended,x_variables);
    % These variables are original decision variables
    x_variables = setdiff(x_variables,aux_variables);
    x_variables = setdiff(x_variables,[extended yalmip('auxvariables')]);

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
    
    
    Graph = yalmip('getdependence');
    Graph = Graph + speye(length(Graph));
    Graph0 = Graph*Graph;
    Graph0 = double(Graph0 | Graph0);
    while ~isequal(Graph,Graph0)    
        Graph = Graph0;
        Graph0 = Graph*Graph;
        Graph0 = double(Graph0 | Graph0);
    end
    
    extended = union(extended, yalmip('auxvariables'));
    % Which auxilliary variables depends only on w?    
    aux_with_only_w_dependence = extended(find(any(Graph(extended,w_variables),2) & ~any(Graph(extended,x_variables),2)));

    % Which auxilliary variables depends on w somehow?    
    aux_with_w_dependence = extended(find((any(Graph(extended,w_variables),2))));    
    aux_with_w_dependence = setdiff(extended(find((any(Graph(extended,w_variables),2)))),aux_with_only_w_dependence);
        
    % Only models x
    aux_without_w_dependence = extended(find(~(any(Graph(extended,w_variables),2))));
        
    x_variables = union(aux_without_w_dependence, x_variables);          
end

allvariables = [depends(F) depends(h)];
x_variables = intersect(allvariables,x_variables);
w_variables = intersect(allvariables,w_variables);

x = recover(x_variables);
w = recover(w_variables);

VariableType.Graph = Graph;
VariableType.x_variables = x_variables;
VariableType.w_variables = w_variables;
VariableType.aux_with_w_dependence = intersect(aux_with_w_dependence,allvariables);
VariableType.aux_with_only_w_dependence = intersect(aux_with_only_w_dependence,allvariables);
