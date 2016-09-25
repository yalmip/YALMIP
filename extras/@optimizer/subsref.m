function varargout = subsref(self,subs)

if isequal(subs.type,'()')
    
    if length(self.diminOrig) > 1
        error('Cannot index in cell-based OPTIMIZER. Perhaps you mean {}');
    end
    if isempty(self.output.z)
        output_length = length(self.output.expression);
    else
        output_length = length(self.map);
    end
    if ~isequal(subs.subs{1},round(subs.subs{1})) || min(subs.subs{1})<1 ||   max(subs.subs{1}) > output_length
        error('Beware of syntax change in optimizer. {} is now used to obtained solution ??? Subscript indices must either be real positive integers or logicals.')
    end
    
    % Create a new function with extracted outputs
    if isempty(self.output.z)
        self.map = self.map(subs.subs{1});
        self.output.expression = self.output.expression(subs.subs{1});
    else
        self.output.expression = self.output.expression(subs.subs{1});
    end
    self.dimoutOrig{1} = size(self.output.expression);
    self.dimout = [numel(self.output.expression) 1];
    varargout{1} = self;
    
elseif isequal(subs.type,'.')
    
    if length(subs) == 1
        switch subs.subs
            case 'options'
                varargout{1} = self.model.options;
            case 'model'
                varargout{1} = self.model;
            otherwise
                error('Field not accesible. You can only acsess P.options')
        end
    else
        error('Field not accesible. You can only acsess P.options')
    end
    
elseif isequal(subs.type,'{}')
    
    if length(subs.subs)>0 && isequal(subs.subs{end},'nosolve')
        NoSolve = 1;
        subs.subs = {subs.subs{1:end-1}};
    else
        NoSolve = 0;
    end
    
    if isa(subs.subs{1},'lmi') || isa(subs.subs{1},'constraint')
        % User instantiates as P{[x == 1, y == ...]}
        cells = cell(1,length(self.diminOrig));
        List = subs.subs{1};
        for i = 1:length(List)
            xi = recover(depends(List(i)));
            for j = 1:length(self.diminOrig)
                if isequal(getvariables(xi), getvariables(self.input.xoriginal{j}))
                    if isequal(getbase(xi),getbase(self.input.xoriginal{j}))
                        B = getbase(List(i));
                        if any(B(:,2:end)>0)
                            cells{j} = full(reshape(-B(:,1),self.diminOrig{j}));
                        else
                            cells{j} = full(reshape(B(:,1),self.diminOrig{j}));
                        end
                    end
                end
            end
        end
        subs.subs = {cells};
    end
    
    if ~isempty(self.ParametricSolution)
        x = subs.subs{1};
        x = x(:);
        [found,j] = isinside(self.ParametricSolution{1}.Pn,x);
        if found
            j = j(1);
            u = self.ParametricSolution{1}.Fi{j}*x + self.ParametricSolution{1}.Gi{j};
            u = reshape(u,self.dimoutOrig{1});
            varargout{1} = u;
            varargout{2} = 0;
        else
            varargout{1} = nan(self.dimoutOrig{1});
            varargout{2} = 1;
        end
        return
    end
    
    if self.model.options.usex0
        if nargout < 5
            warning('If you intend to use initial guesses, you must save fifth output [sol,~,~,~,P] = P{p}');
        end
    end
    
    % This is not really supported yet...
    if isa(subs.subs{1},'sdpvar')
        varargout{1} = yalmip('definemulti','optimizer_operator',subs(1).subs{1},self,self.dimout);
        return
    end
    
    if ~isempty(self.diminOrig)
        % Normalize to a cell-format
        if ~isa(subs.subs{1},'cell')
            subs.subs{1} = {subs.subs{1}};
        end
        
        % Check number of cells
        if length(subs.subs{1})~=length(self.diminOrig) && ~isempty(self.diminOrig)
            error('The number of cell elements in input does not match OPTIMIZER declaration');
        end
        
        % Realify...
        for i = 1:length(subs.subs{1})
            if self.complexInput(i)
                subs.subs{1}{i} = [real(subs.subs{1}{i});imag(subs.subs{1}{i})];
            end
        end
        
        
        % If blocked call, check so that there are as many blocks for every
        % argument. Note, we only analyze instantiated blocks
        dimBlocks = nan(length(subs.subs{1}),1);
        for i = 1:length(subs.subs{1})
            dimBlocks(i) = numel(subs.subs{1}{i}) / prod(self.diminOrig{i});
        end
        if any(dimBlocks)>1 && any(dimBlocks==0)
            error('Blocked data in partial instantiation is not possible');
        end
        if all(dimBlocks)
            if ~isnan(dimBlocks)
                if ~all(dimBlocks == dimBlocks(i))
                    error('Dimension mismatch on the input arguments compared to definition');
                end
                if any(dimBlocks ~= fix(dimBlocks))
                    error('Dimension mismatch on the input arguments compared to definition');
                end
            end
        end
  
        % Just pick out any element, it the number of blocks, as all should
        % be the same
        nBlocks = max(dimBlocks(find(dimBlocks)));
                
        left = ones(1,length(subs.subs{1}));
        aux = [];
        aux2 = [];
        suppliedData = [];
        for i = 1:nBlocks
            aux2 = [];
            try
                for j = 1:length(subs.subs{1})
                    data = subs.subs{1}{j};
                    if isempty(data)
                        temp = nan(self.diminOrig{j});
                        temp = temp(:);
                        suppliedData(j) = 0;
                    else
                        temp = data(:,left(j):self.diminOrig{j}(2)+left(j)-1,:);
                        suppliedData(j) = 1;
                    end
                    temp = temp(:);
                    temp = temp(self.mask{j});
                    aux2 = [aux2;temp];
                    left(j) = left(j) + self.diminOrig{j}(2);
                end
            catch
                error('The dimension on a supplied parameter value does not match the definition.');
            end
            %left = left + self.diminOrig{1}(2);
            aux = [aux aux2];
        end
        subs.subs{1} = aux;
        
        % Input is now given as [x1 x2 ... xn]
        u = [];
      %  nBlocks = dimBlocks(1);
        start = 1;
        if isnan(nBlocks)
            nBlocks = 1;
        end
    else
        u = [];
        start = 0;
        nBlocks = 1;
    end
    
    for i = 1:nBlocks
        if ~isempty(self.diminOrig)
            thisData = subs.subs{1}(:,start:start + self.dimin(2)-1);
        else
            thisData = [];
        end
        if NoSolve || any(isnan(thisData)) || any(self.instatiatedvalues) || self.nonlinear && ~self.complicatedEvalMap%isempty(self.model.evalMap)
            originalModel = self.model;
             
            if any(isnan(thisData))
                % Partial elimination)
                current_parametric = self.parameters(~isnan(thisData));
                all_parametric = self.parameters;
                % Crate an object with a reduced set of variables to
                % eliminate
                self.parameters = current_parametric;
                self = optimizer_precalc(self);
                % Eliminate
                thisData = thisData(~isnan(thisData));
                [self.model,keptvariables,infeasible] = eliminatevariables(self.model, current_parametric,thisData(:),all_parametric);
                % Update as new optimizer in remaining variables, remap to
                % currently used variables
                old_global_parametric_index = find(ismember(self.orginal_usedvariables,self.model.used_variables(current_parametric)));
                self.instatiatedvalues(old_global_parametric_index) = thisData;
                self.parameters = find(ismember(keptvariables,setdiff(all_parametric,current_parametric)));
                self.model.used_variables = self.model.used_variables(keptvariables);
                self.dimin = [length(self.parameters) 1];
                self.diminOrig = {self.diminOrig{find(~suppliedData)}};
                self.input.xoriginal = {self.input.xoriginal{find(~suppliedData)}};
                self.mask = {self.mask{find(~suppliedData)}};
                self = optimizer_precalc(self);
                varargout{1} = self;
                
                return
            else
                % Standard case where we will solve a problem
                [self.model,keptvariables,infeasible] = eliminatevariables(self.model,self.parameters,thisData(:),self.parameters);
            end
            
            if isnan(infeasible)
                % Elimination was not run as there were no variables
                % Data is inherited from last eliminate
                infeasible = self.infeasible;
                keptvariables = self.keptvariables;
            end
            % Turn off equality presolving for simple programs. equality
            % presolve has benefits when the are stuff like log
            self.model.presolveequalities = length(self.model.evalMap > 0);
            if ~infeasible
                if self.model.options.usex0 && ~isempty(self.lastsolution)
                    self.model.x0 = zeros(length(self.model.c),1);
                    self.model.x0 = self.lastsolution;
                elseif ~self.model.options.usex0
                    self.model.x0 = [];
                end
                if NoSolve
                    % We just instantiate the model, and return it
                    self.dimin = [];
                    self.parameters = [];
                    self.diminOrig = {};
                    self.infeasible = infeasible;
                    self.keptvariables = keptvariables;
                    varargout{1} = self;
                    return
                else
                    eval(['output = ' self.model.solver.call '(self.model);']);
                end
                if output.problem == 0 && self.model.options.usex0
                    self.lastsolution = output.Primal;
                end
                x = self.instatiatedvalues;             
                x(ismember(self.orginal_usedvariables,self.model.used_variables(self.parameters))) = thisData(:);
                x(find(ismember(self.orginal_usedvariables,self.model.used_variables(keptvariables)))) = output.Primal;
                %x = zeros(length(self.orginal_usedvariables),1);%originalModel.c*0;                
                %x(self.parameters) = thisData(:); % Bugfix aus (https://groups.google.com/forum/#!category-topic/yalmip/S0ukzE_wLGs)              
                %x(keptvariables) = output.Primal;
                output.Primal = x;
            else
                output.problem = 1;
                output.Primal = originalModel.c*0;
                output.Dual = [];
            end
            originalModel.precalc = self.model.precalc;
            self.model = originalModel;
        else
            if ~isempty(thisData)
                if ~isempty(self.model.evalMap) &&  ~self.model.solver.evaluation
                    error('After fixing parameters, there are still nonlinear operators in the model, but the solver does not support this. Note that YALMIP is not guaranteed to remove all operators, even though they only contain parametric expressions. As an example, exp(1+parameter) will not be reduced, while exp(parameter) will. You will have to use another solver, or reparameterize your model (look at exp(1+parameter) as a new parameter instead)');
                end
                self.model.F_struc(1:prod(self.dimin),1) = thisData(:);
            end
            if self.model.options.usex0 && ~isempty(self.lastsolution)
                self.model.x0 = self.lastsolution;
            elseif ~self.model.options.usex0
                self.model.x0 = [];
            end
            
            if NoSolve
                % We just instantiate the model, and return it
                self.dimin = [];
                self.parameters = [];
                self.diminOrig = {};
                self.infeasible = 0;
                self.keptvariables = 1:length(self.model.c);
                varargout{1} = self;
                return
            else
                output = self.model.solver.callhandle(self.model);
            end
            if output.problem == 0 && self.model.options.usex0
                self.lastsolution = output.Primal;
            end
        end
        if output.problem==1
            output.Primal = output.Primal+nan;
        end
        if isempty(self.output.z)
            if ~isempty(output.Primal)
                % Make sure we map index 0 to Nans
                % 0 corresponds to variables which weren't visible in
                % problem
                output.Primal = [nan;output.Primal];
                u = [u reshape(output.Primal(1+self.map),self.dimout)];
            else
                u = [u reshape(0*self.map+nan,self.dimout)];
            end
        else
            if ~isempty(output.Primal)
                % Make sure we map index 0 to Nans
                % 0 corresponds to variables which weren't visible in
                % problem
                output.Primal = [nan;output.Primal];
                assign(self.output.z,output.Primal(1+self.map));
                assign(self.input.expression,thisData);
                u = [u reshape(double(self.output.expression),self.dimout)];
            end
        end
        varargout{2}(i) = output.problem;
        varargout{3}{i} = yalmiperror(output.problem);
        varargout{4}{i} = output.Dual;
        if ~isempty(self.dimin)
            start = start + self.dimin(2);
        end
    end
    if length(self.dimoutOrig)>1
        % top = 1;
        realDimOut = self.dimoutOrig;
        allu = cell(1, length(self.dimoutOrig));
        for k = 1:nBlocks
            top = 1;
            for i = 1:length(self.dimoutOrig)
                n = prod(self.dimoutOrig{i});
                uvec = reshape(u(top:top+n-1,k),self.dimoutOrig{i});
                if self.complexOutput(i)
                    uvec = uvec(1:size(uvec,1)/2,:) +  uvec(1+size(uvec,1)/2:end,:)*sqrt(-1);
                end
                allu{i} = [allu{i} uvec];
                top = top + n;
            end
        end
        varargout{1} = allu;
    elseif nBlocks==1
        varargout{1} = reshape(u(:),self.dimoutOrig{1});
        if self.complexOutput(1)
            if length(self.dimoutOrig{1})==2
                varargout{1} = varargout{1}(1:self.dimoutOrig{1}(1)/2,:) + sqrt(-1)*varargout{1}(self.dimoutOrig{1}(1)/2+1:end,:);
            else
                % FIXME: This only works for 3D...
                varargout{1} = varargout{1}(1:self.dimoutOrig{1}(1)/2,:,:) + sqrt(-1)*varargout{1}(self.dimoutOrig{1}(1)/2+1:end,:,:);
            end
        end
    else
        varargout{1} = u;
    end
    varargout{5} = self;
    varargout{6} = output;
end
