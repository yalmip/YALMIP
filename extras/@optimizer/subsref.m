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
    if ~isequal(subs.subs{1},round(subs.subs{1})) | min(subs.subs{1})<1 |   max(subs.subs{1}) > output_length
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
    self.dimout = [prod(size(self.output.expression)) 1];
    varargout{1} = self;
    
elseif isequal(subs.type,'.')
    
    switch subs.subs
        case 'model'
            varargout{1} = self.model;
        otherwise
            error('Field not accesible')
    end
    
elseif isequal(subs.type,'{}')

    % This is not really supported yet...
    if isa(subs.subs{1},'sdpvar')
        varargout{1} = yalmip('definemulti','optimizer_operator',subs(1).subs{1},self,self.dimout);
        return
    end
    
    % Normalize to a cell-format
    if ~isa(subs.subs{1},'cell')
        subs.subs{1} = {subs.subs{1}};
    end
    
    % Check number of cells
    if length(subs.subs{1})~=length(self.diminOrig)
        error('The number of cell elements in input does not match OPTIMIZER declaration');
    end
          
    % Realify...
    for i = 1:length(subs.subs{1})
        if self.complexInput(i)
            subs.subs{1}{i} = [real(subs.subs{1}{i});imag(subs.subs{1}{i})];
        end
    end
    
    % If blocked call, check so that there are as many blocks for every
    % argument
    for i = 1:length(subs.subs{1})
        dimBlocks(i) = prod(size(subs.subs{1}{i})) / prod(self.diminOrig{i});
    end    
    if ~isnan(dimBlocks) & ~all(dimBlocks == dimBlocks(i))
        error('Dimension mismatch on the input arguments compared to definition');
    end
    if ~isnan(dimBlocks)
        if any(dimBlocks ~= fix(dimBlocks))
            error('Dimension mismatch on the input arguments compared to definition');
        end
    end
    
    
    left = ones(1,length(subs.subs{1}));
    aux = [];
    aux2 = [];
    for i = 1:dimBlocks(1)
        aux2 = [];
        for j = 1:length(subs.subs{1})
            temp = subs.subs{1}{j}(:,left(j):self.diminOrig{j}(2)+left(j)-1,:);
            temp = temp(:);
            temp = temp(self.mask{j});
            aux2 = [aux2;temp];
             left(j) = left(j) + self.diminOrig{j}(2);
        end
        %left = left + self.diminOrig{1}(2);
        aux = [aux aux2];
    end
    subs.subs{1} = aux;
    
    % Input is now given as [x1 x2 ... xn]   
    u = [];
    nBlocks = dimBlocks(1);
    start = 1;
    if isnan(nBlocks)
        nBlocks = 1;
    end
        
    for i = 1:nBlocks
        thisData = subs.subs{1}(:,start:start + self.dimin(2)-1);
        if self.nonlinear & ~self.complicatedEvalMap%isempty(self.model.evalMap)
            originalModel = self.model;
            [self.model,keptvariables,infeasible] = eliminatevariables(self.model,self.parameters,thisData(:));
            % Turn off equality presolving for simple programs. equality
            % presolve has benefits when the are stuff like log
            self.model.presolveequalities = length(self.model.evalMap > 0);
            if ~infeasible                          
                eval(['output = ' self.model.solver.call '(self.model);']);
                x = originalModel.c*0;
                x(keptvariables) = output.Primal;
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
                if ~isempty(self.model.evalMap) &  ~self.model.solver.evaluation
                    error('After fixing parameters, there are still nonlinear operators in the model, but the solver does not support this. Note that YALMIP is not guaranteed to remove all operators, even though they only contain parametric expressions. As an example, exp(1+parameter) will not be reduced, while exp(parameter) will. You will have to use another solver, or reparameterize your model (look at exp(1+parameter) as a new parameter instead)');
                end
                self.model.F_struc(1:prod(self.dimin),1) = thisData(:);
            end         
            output = self.model.solver.callhandle(self.model);           
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
        start = start + self.dimin(2);
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
end