function varargout = subsref(self,subs)

if isequal(subs.type,'()')
    
    % New syntax. () used to be replacement/computation
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
        self.map = self.map(subs.subs{1});
        self.output.expression = self.output.expression(subs.subs{1});
    end
    self.dimout = size(self.output.expression);
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
        u = yalmip('addextendedvariable','optimizer_operator',self,subs.subs{1});
        return
    end

    % Input is given as [x1 x2 ... xn]
    if size(subs.subs{1},1) == self.dimin(1)
        % Check that width is an multiple of parameter width
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
    for i = 1:nBlocks
        thisData = subs.subs{1}(:,start:start + self.dimin(2)-1);
        %        if self.nonlinear & isempty(self.model.evalMap) & isempty(self.model.bilinear_variables) & isempty(self.model.integer_variables)
        if self.nonlinear & isempty(self.model.evalMap) & isempty(self.model.integer_variables)
            originalModel = self.model;
            [self.model,keptvariables,infeasible] = eliminatevariables(self.model,self.parameters,thisData(:));
            if ~infeasible
                eval(['output = ' self.model.solver.call '(self.model);']);
                x = originalModel.c*0;
                x(keptvariables) = output.Primal;
                output.Primal = x;
            else
                output.problem = 1;
                output.Primal = originalModel.c*0;
            end
        else
            
            %[ii,jj,kk] = find(self.model.F_struc(1:prod(self.dimin),2:end))
            
            if ~isempty(thisData)
                self.model.F_struc(1:prod(self.dimin),1) = thisData(:);
            end         
            output = self.model.solver.callhandle(self.model);
           % eval(['output = ' self.model.solver.call '(self.model);']);
        end
        if output.problem==1
            output.Primal = output.Primal+nan;
        end
        if isempty(self.output.z)
            if ~isempty(output.Primal)
                u = [u output.Primal(self.map)];
            else
                u = [u 0*self.map+nan];
            end
        else
            if ~isempty(output.Primal)
                assign(self.output.z,output.Primal(self.map));
                u = [u double(self.output.expression)];
            end
        end
        varargout{2}(i) = output.problem;
        varargout{3}{i} = yalmiperror(output.problem);
        varargout{4}{i} = output.Dual;
        start = start + self.dimin(2);
    end
    varargout{1} = u;
end