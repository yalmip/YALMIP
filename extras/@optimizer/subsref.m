function varargout = subsref(self,subs)

if isequal(subs.type,'()')
    
    % Previously not allowed and user was to use cell format for some now
    % forgotten reason (ideas about other functionalities for () perhaps)
    % We support it now though, so throw to the previously required {}
    subs.type ='{}';
    [varargout{1:nargout}] = subsref(self,subs);
    
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
       
    % User has typed X();
    if isempty(subs.subs)
        subs.subs = {[]};
    else
        % Flatten the expression to account for various cell calls
        aux = {};
        for i = 1:length(subs.subs)
            if isa(subs.subs{i},'cell')
                for j = 1:length(subs.subs{i});
                    aux = {aux{:},subs.subs{i}{j}};
                end
            else
                aux = {aux{:},subs.subs{i}};
            end
        end
        subs.subs = aux;
    end
    
    if length(self.dimin)==0
        if length(subs.subs)==1
            if ~isempty(subs.subs{1})
                error('No parameters left to assign.');
            end
        elseif length(subs.subs)>1
            error('No parameters left to assign.');
        end
    end
    
    NoSolve = 0;
    % Check for nosolve flag
    if length(subs.subs)>0 && isequal(subs.subs{end},'nosolve')
        NoSolve = 1;
        subs.subs = {subs.subs{1:end-1}};   
    end
          
    if isa(subs.subs{1},'lmi') || isa(subs.subs{1},'constraint')
        % User instantiates as P{[x == 1, y == ...]}
        cells = cell(1,length(self.diminOrig));
        for i = 1:length(cells)
            cells{i} = []; % Default user has not specified, partial inst.
        end
        List = subs.subs{1};
        for i = 2:length(subs.subs);List = [List, subs.subs{i}];end        
        matched = zeros(1,length(List));
        for i = 1:length(List)
            xi = recover(depends(List(i)));
            for j = 1:length(self.diminOrig)
                if isequal(getvariables(xi), getvariables(self.input.xoriginal{j}))
                    if isequal(getbase(xi),getbase(self.input.xoriginal{j}))
                        B = getbase(List(i));
                        if any(B(:,2:end)>0)
                            matched(i) = 1;
                            cells{j} = full(reshape(-B(:,1),self.diminOrig{j}));
                        else
                            matched(i) = 1;
                            cells{j} = full(reshape(B(:,1),self.diminOrig{j}));
                        end
                    else
                        error('Only simple expressions allowed when declaring variable values')
                    end
                end
            end
        end
        if ~all(matched)
            error('Only simple expressions allowed. Could not match assignments in optimizer call');
        end
        subs.subs = cells;
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
    
    if self.model.options.warmstart
        if nargout < 5
            warning('If you intend to use initial guesses, you must use a fifth output as [sol,problem,~,~,P] = P(p)');
        end
    end
    
    % This is not really supported yet...
    if isa(subs.subs{1},'sdpvar')
        varargout{1} = yalmip('definemulti','optimizer_operator',subs(1).subs{1},self,self.dimout);
        return
    end
    
    if ~isempty(self.diminOrig)
       
        % Check number of cells
        if length(subs.subs)~=length(self.diminOrig) && ~isempty(self.diminOrig)
            error('The number of cell elements in input does not match OPTIMIZER declaration');
        end
        
        % Realify...
        for i = 1:length(subs.subs)
            if self.complexInput(i)
                subs.subs{i} = [real(subs.subs{i});imag(subs.subs{i})];
            end
        end
        
        % Check for garbage call
        for i = 1:length(subs.subs)
            if isnan(subs.subs{i})
                if length(self.dimoutOrig)>1
                    for j = 1:length(self.dimoutOrig)
                        sol{j} = nan(self.dimoutOrig{j});
                    end
                else
                    sol = nan(self.dimoutOrig{1});
                end
                varargout{1} = sol;
                varargout{2} = -10;
                varargout{3} = yalmiperror(-10);
                varargout{4} = []; 
                varargout{5} = self; 
                varargout{6} = self; 
                return
            end
        end
        
        
        % If blocked call, check so that there are as many blocks for every
        % argument. Note, we only analyze instantiated blocks
        dimBlocks = nan(length(subs.subs),1);
        for i = 1:length(subs.subs)
            dimBlocks(i) = numel(subs.subs{i}) / prod(self.diminOrig{i});
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
        
        left = ones(1,length(subs.subs));
        aux = [];
        aux2 = [];
        suppliedData = [];
        for i = 1:nBlocks
            aux2 = [];
            try
                for j = 1:length(subs.subs)
                    data = subs.subs{j};
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
        subs.subs = aux;
        
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
            thisData = subs.subs(:,start:start + self.dimin(2)-1);
        else
            thisData = [];
        end
        
        originalModel = self.model;
        
        if any(isnan(thisData))
            % Partial elimination)
            currentParametricIndex = self.model.parameterIndex(~isnan(thisData));
            thisData = thisData(~isnan(thisData));
                        
            allParametricIndex = self.model.parameterIndex;
            globalParametricIndex = find(ismember(self.orginal_usedvariables,self.model.used_variables(currentParametricIndex)));
            self.instatiatedvalues(globalParametricIndex) = thisData;
            
            evalParameters = [];
            for k = 1:length(self.model.evalMap)
                for j = 1:length(currentParametricIndex)
                    if isequal(self.model.evalMap{k}.variableIndex,currentParametricIndex(j))
                        evalParameters = [evalParameters self.model.evalMap{k}.computes];
                    end
                end
            end
            self.model.evalParameters = evalParameters;
            
            
            % Crate an object with a reduced set of variables to
            % eliminate
            self.model.parameterIndex = currentParametricIndex;
            self = optimizer_precalc(self);                 
            [self.model,keptvariablesIndex] = eliminatevariables(self.model,currentParametricIndex,thisData(:),allParametricIndex);
            
            % Update as new optimizer in remaining variables, remap
            % parameters to currently used variables                                   
            remainingParametersIndex = setdiff(allParametricIndex,currentParametricIndex,'stable');            
            [~,loc] = ismember(remainingParametersIndex,keptvariablesIndex);
            self.model.parameterIndex = loc;
            
            % Remap all evaluation operators
            for k = 1:length(self.model.evalMap)
                self.model.evalMap{k}.computes = find(self.model.evalMap{k}.computes == keptvariablesIndex);
                self.model.evalMap{k}.variableIndex = find(self.model.evalMap{k}.variableIndex == keptvariablesIndex);
            end
            
            evalParameters = [];
            for k = 1:length(self.model.evalMap)
                for j = 1:length(self.model.parameterIndex)
                    if isequal(self.model.evalMap{k}.variableIndex,self.model.parameterIndex(j))
                        evalParameters = [evalParameters self.model.evalMap{k}.computes];
                    end
                end
            end
            self.model.evalParameters = evalParameters;
                                    
            self.dimin = [length(self.model.parameterIndex) 1];                     
            self.diminOrig = {self.diminOrig{find(~suppliedData)}};
            self.input.xoriginal = {self.input.xoriginal{find(~suppliedData)}};
            self.mask = {self.mask{find(~suppliedData)}};
            self = optimizer_precalc(self);
            varargout{1} = self;
            return
        else                                   
            % Standard case where we eliminate all variables left
            self.instatiatedvalues(ismember(self.orginal_usedvariables,self.model.used_variables(self.model.parameterIndex))) = thisData(:);
            [self.model,keptvariablesIndex] = eliminatevariables(self.model,self.model.parameterIndex,thisData(:),self.model.parameterIndex);
            
            % Remap all evaluation operators
            self.model.evalVariables = [];
            for k = 1:length(self.model.evalMap)
                self.model.evalMap{k}.computes = find(self.model.evalMap{k}.computes == keptvariablesIndex);
                self.model.evalMap{k}.variableIndex = find(self.model.evalMap{k}.variableIndex == keptvariablesIndex);
                self.model.evalVariables = [self.model.evalVariables self.model.evalMap{k}.computes];
            end                                         
            self.model.evalVariables = sort(self.model.evalVariables);            
        end

        % Turn off equality presolving for simple programs. equality
        % presolve has benefits when the are stuff like log
        self.model.presolveequalities = length(self.model.evalMap) > 0;
        if ~self.model.infeasible
            if self.model.options.warmstart && ~isempty(self.lastsolution)
                self.model.x0 = zeros(length(self.model.c),1);
                self.model.x0 = self.lastsolution;
            elseif ~self.model.options.warmstart
                self.model.x0 = [];
            end
            if NoSolve
                % We just returns the presolved model
                self.model.parameterIndex = [];
                self.dimin = [];                
                self.diminOrig = {};                         
                varargout{1} = self;
                return
            else
                eval(['output = ' self.model.solver.call '(self.model);']);
            end
            
            if output.problem == 0 && self.model.options.warmstart
                self.lastsolution = output.Primal;
            end
            x = self.instatiatedvalues;
            if isempty(output.Primal)
                v = find(ismember(self.orginal_usedvariables,self.model.used_variables));
                x(v) = nan(length(v),1);
            else
                x(ismember(self.orginal_usedvariables,self.model.used_variables)) = output.Primal;
            end
            output.Primal = x;
            
        else
            output.problem = 1;
            output.Primal = originalModel.c*0;
            output.Dual = [];
        end
        originalModel.precalc = self.model.precalc;
        self.model = originalModel;
        
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
                if ~isempty(thisData)
                    assign(self.input.expression,thisData);
                end
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
