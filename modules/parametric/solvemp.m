function varargout = solvemp(F,h,ops,x,y)
%SOLVEMP Computes solution to multi-parametric optimization problem
%
% min_z(x)   h(x,z)
% subject to
%            F(x,z) >= 0
%
%
% [SOL, DIAGNOSTIC,Z,HPWF,ZPWF] = SOLVEMP(F,h,options,x,y)
%
% SOL        : Multi-parametric solution (see MPT toolbox)
%
% DIAGNOSTIC : struct with diagnostic information
%
% Z          : SDPVAR object with the detected decision variable z
%
% HPWF       : The value function as a pwf function
%
% ZPWF       : The optimal decision variable as a pfw function
%
% Input
%    F        : Object describing the constraints.
%    h        : SDPVAR object describing the objective function h(x,z).
%    options  : solver options. See SDPSETTINGS. Can be [].
%    x        : Parametric variables
%    y        : Requested decision variables (subset of z)
%
% NOTE : If you are solving a problem leading to an mpMILP, the
% output SOL will be a set-valued map. To obtain the minimal
% solution (without so called overlaps), run removeOverlaps(SOL). If you
% have requested the 5th output ZPWF, overlaps are automatically removed.
% If your problem leads to an mpMIQP,  the output SOL will also be a
% set-valued map, but there is currently no way in MPT to obtain a
% non-overlapping solution. To use the solution in MPT, the command
% mpt_mergeCS(SOL) can be useful. Notice that the fifth output argument
% not will be available for mpMIQP problems.
%
% See also PARAMETRIC, SET, SDPSETTINGS, YALMIPERROR

if nargin <= 3
    ops = sdpsettings;
end

if nargin <=3
    x = [];
    y = [];
end

if isa(F,'constraint')
    F = lmi(F);
end

par_declarations = is(F,'parametric');
if any(par_declarations)
    x = [x;recover(getvariables(sdpvar(F(find(par_declarations)))))];
    F = F(find(~par_declarations));
end

if length(x) == 0
    error('solvemp must always have 4 input arguments or a parametric declaration');
end

if ~isempty(ops)
    if isequal(ops.solver,'')
        ops.solver = 'mpt,pop';
    end
else
    ops = sdpsettings('solver','mpt,pop');
end

if nargin == 4
    y = [];
    ny = 0;
    my = 0;
else
    % YALMIP wants a vector as desired decsision variable
    [ny,my] = size(y);
    y = reshape(y,ny*my,1);
end

% Robustify first?
if length(F) > 0
    unc_declarations = is(F,'uncertain');
    if any(unc_declarations)
        w = recover(getvariables(sdpvar(F(find(unc_declarations)))));
        F = F(find(~unc_declarations));
        [F,h,failure] = robustify(F,h,ops,w);
        if failure
            error('Derivation of robust counter-part failed')
        end
    end
end

if max(size(h))>1
    error('Objective function must be scalar or empty');
end

sol = solvesdp(F,h,ops,x,y);

if isfield(sol,'mpsol')
    if ~isfield(sol.mpsol,'model')
        varargout{1} = [];
        varargout{2} = sol;
        varargout{3} = [];
        varargout{4} = [];
        varargout{5} = [];
    elseif isempty(sol.mpsol.model{1})
        varargout{1} = sol.mpsol.model;       
        varargout{2} = sol;
        varargout{3} = [];
        varargout{4} = [];
        varargout{5} = []; 
    else

        mpsolution = sol.mpsol.model;
        varargout{1} = sol.mpsol.model;

        if nargout > 2
            z = recover(sol.solveroutput.U);
            x = recover(sol.solveroutput.x);
            varargout{3}= z;
        end

        if nargout > 3
            % User wants the value function                
            if length(mpsolution) == 1
                if isequal(mpsolution{1}.convex,1)
                    % Simple mpLP value function
                    if ops.mp.simplify
                        s = mpsolution{1};
                        s.Fi = s.Bi;
                        s.Gi = s.Ci;
                        s = mpt_simplify(s);
                        s.Bi = s.Fi;
                        s.Ci = s.Gi;
                        varargout{4} = pwf(s,x,'convex');
                    else
                        varargout{4} = pwf(mpsolution{1},x,'convex');
                    end
                else
                    % Probably generated from removing overlaps
                    varargout{4} = pwf(mpsolution,x,'general');
                end
            else
                % No overlap removal done
                varargout{4} = pwf(mpsolution,x,'convexoverlapping');
            end
        end

        if nargout > 4
            % User wants optimizer in YALMIP format
            % Any overlaps?
            anyQP = 0;
            if length(varargout{1}) > 1
                for i = 1:length(sol.mpsol.model)
                    if nnz([sol.mpsol.model{i}.Ai{:}])>0
                        anyQP = 1;
                        break
                    end
                end
                if ~anyQP
                    minimalmodel{1} = mpt_removeOverlaps(sol.mpsol.model);
                    varargout{1} = minimalmodel;
                end
            else
                minimalmodel = varargout{1};
            end
            % PWA assumes we want Bi and Ci
            if ~anyQP
                minimalmodel{1}.valuefunction.Bi =  minimalmodel{1}.Bi;
                minimalmodel{1}.valuefunction.Ci =  minimalmodel{1}.Ci;
                minimalmodel{1}.Ai = cell(1,length(minimalmodel{1}.Fi));
                minimalmodel{1}.Bi = minimalmodel{1}.Fi;
                minimalmodel{1}.Ci = minimalmodel{1}.Gi;
                varargout{5} = pwf(minimalmodel,x,'general');
                if min([ny my])>0
                    varargout{5} = reshape(varargout{5},ny,my);
                end
            else
                disp('Optimizer (5th output) not available for overlapping quadratic problems.');
                varargout{5} = [];
            end
        end
    end
else
    varargout{1} = [];
    varargout{2} = sol;
    varargout{3} = [];
    varargout{4} = [];
    varargout{5} = [];
end
varargout{2} = sol;

