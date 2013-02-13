function varargout = subsref(self,subs)

if isequal(subs.type,'()')
    
    if size(self.diminOrig,1) > 1
        error('Can index in cell-based OPTIMIZER. Perhaps you mean {}');
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
     %   self.map = self.map(subs.subs{1});
        self.output.expression = self.output.expression(subs.subs{1});
    end
    self.dimoutOrig = size(self.output.expression);
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
        
    if size(self.diminOrig,1) > 1
        % Concatenate inputs?
        if isa(subs.subs{1},'cell')
            if length(subs.subs{1})~=size(self.diminOrig,1)
                error('The number of cell elements in input does not match OPTIMIZER declaration');
            end
            vecIn = [];
            for i = 1:length(subs.subs{1})
                arg = subs.subs{1}{i};
                if ~isequal(size(arg),self.diminOrig(i,:))
                    error('Size mismatch on argument in cell');
                end
                temp = arg(:);
                temp = temp(self.mask{i});
                vecIn = [vecIn;temp];
            end    
            subs.subs{1} = vecIn;
        else
            error('OPTIMIZER defined with cell input format, but you sent a vector');
        end
    else
        if isa(subs.subs{1},'cell')
            if length(subs.subs{1})>1
                error('You are sending multiple cells, but OPTIMIZER was defined with only one argument');
            else
                subs.subs{1} = subs.subs{1}{1};
            end
        end 
        if ~isequal(size(subs.subs{1},1),self.diminOrig(1))
            error('Size mismatch. The parameter does match the size used in definition of OPTIMIZER object')
        end
        subs.subs{1} =  subs.subs{1}(:);
        subs.subs{1} = reshape(subs.subs{1},prod(self.diminOrig),[]);
        subs.subs{1} = subs.subs{1}(self.mask{1},:);
    end

    % Input is now given as [x1 x2 ... xn]
    % Note: this is not supported in cell based arguments
    % If cell-based argument was used, only one call can be made at once
    % i.e., n=1
    if size(subs.subs{1},1) == self.dimin(1)
        % Check that wP2idth is an multiple of parameter width
        if mod(size(subs.subs{1},2),self.dimin(2))>0
            error('Input argument has wrong size (The width is not a multiple of parameter width');
        end
    else
        error('Input argument has wrong size (The height does not match with you original argument)');
    end

    u = [];
    nBlocks = size(subs.subs{1},2)/self.dimin(2);
    start = 1;
    if isnan(nBlocks)
        nBlocks = 1;
    end
    
    if nBlocks>1 && (size(self.dimoutOrig,1)>1 || size(self.diminOrig,1)>1)
        error('Multiple calls not supported in OPTIMIZER when using cell format');
    end
    
    for i = 1:nBlocks
        thisData = subs.subs{1}(:,start:start + self.dimin(2)-1);
        if self.nonlinear & isempty(self.model.evalMap)
            originalModel = self.model;
            try
                [self.model,keptvariables,infeasible] = eliminatevariables(self.model,self.parameters,thisData(:));
            catch
                error('Nonlinear replacement in optimizer object only supported in MATLAB R2012A or later');
            end
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
            originalModel.precalc.S = self.model.precalc.S;
            originalModel.precalc.skipped = self.model.precalc.skipped;
            originalModel.precalc.newmonomtable = self.model.precalc.newmonomtable;
            self.model = originalModel;
        else                                  
            if ~isempty(thisData)
                self.model.F_struc(1:prod(self.dimin),1) = thisData(:);
            end         
            output = self.model.solver.callhandle(self.model);           
        end
        if output.problem==1
            output.Primal = output.Primal+nan;
        end
        if isempty(self.output.z)
            if ~isempty(output.Primal)
                u = [u reshape(output.Primal(self.map),self.dimout)];
            else
                u = [u reshape(0*self.map+nan,self.dimout)];
            end
        else
            if ~isempty(output.Primal)
                assign(self.output.z,output.Primal(self.map));
                u = [u reshape(double(self.output.expression),self.dimout)];
            end
        end
        varargout{2}(i) = output.problem;
        varargout{3}{i} = yalmiperror(output.problem);
        varargout{4}{i} = output.Dual;       
        start = start + self.dimin(2);
    end
    if size(self.dimoutOrig,1)>1    
        top = 1;
        realDimOut = self.dimoutOrig;
        allu = {};
        for i = 1:size(realDimOut,1)
            allu{i} = reshape(u(top:top+realDimOut(i,1)*realDimOut(i,2)-1),realDimOut(i,1),realDimOut(i,2));
            top = top + realDimOut(i,1)*realDimOut(i,2);
        end
        varargout{1} = allu;
    elseif nBlocks==1
        varargout{1} = reshape(u(:),self.dimoutOrig);
    else
        varargout{1} = u;
    end  
    varargout{5} = self;
end