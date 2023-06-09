function options = selectSOSmodel(F,options,NonLinearParameterization,noRANK,IntegerData,UncertainData,obj)
switch options.sos.model
    case 0
        constraint_classes = constraintclass(F);
        noCOMPLICATING = ~any(ismember([7 8 9 10 12 13 14 15],constraint_classes));
        if noCOMPLICATING & ~NonLinearParameterization && noRANK && ~IntegerData && ~isa(obj,'logdet')
            options.sos.model = 1;
            if options.verbose>0;disp('Using kernel representation (options.sos.model=1).');end
        else
            if NonLinearParameterization
                if options.verbose>0;disp('Using image representation (options.sos.model=2). Nonlinear parameterization found');end
            elseif ~noRANK
                if options.verbose>0;disp('Using image representation (options.sos.model=2). SOS-rank constraint was found.');end
            elseif IntegerData
                if options.verbose>0;disp('Using image representation (options.sos.model=2). Integrality constraint was found.');end
            elseif UncertainData
                if options.verbose>0;disp('Using image representation (options.sos.model=2). Uncertain data was found.');end                
            elseif isa(obj,'logdet')
                if options.verbose>0;disp('Using image representation (options.sos.model=2). Logdet objective was found.');end 
            else
                if options.verbose>0;disp('Using image representation (options.sos.model=2). Integer data, KYPs or similar was found.');end
            end
            options.sos.model = 2;
        end
    case 1
        if NonLinearParameterization
            if options.verbose>0;disp('Switching to image model due to nonlinear parameterization (not supported in kernel model).');end
            options.sos.model = 2;
        end
        if ~noRANK
            if options.verbose>0;disp('Switching to image model due to SOS-rank constraints (not supported in kernel model).');end
            options.sos.model = 2;
        end
        if IntegerData
            if options.verbose>0;disp('Switching to image model due to integrality constraints (not supported in kernel model).');end
            options.sos.model = 2;
        end        
    case 3
    otherwise
end
