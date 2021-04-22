function p = lastresortBounds(p,upper)

p.boundGenerator.optimizer = [];
if p.options.bmibnb.sdpbounder && isempty(p.evalMap) && all(p.variabletype <=2)
    % Only do this if very badly bounded, or user really wants it
    if (all(isinf(p.ub)) || all(isinf(p.lb))) || p.options.bmibnb.sdpbounder>=1
        if ~isempty(p.solver.sdpsolver)
            % QP with all unbounded variables. This never ends well.
            % Try to solve a semidefinite relaxation to obtain an upper
            % bound on x'*x. 
            % Currently only implemented when there are linear constraints,
            % can be generalized in absurdum. x'*x <= gamma when
            % f+c'*x+x'*Q*x <= U, Ax == b  
            % [gamma 0;0 I] - tau*[U-f -c/2;-c/2 -Q] - tau2*[b;-A]'*[b -A] >= 0
            
            b = p.F_struc(1:p.K.f,1);
            A = p.F_struc(1:p.K.f,2:end);
            v = find(p.variabletype == 0);
            dontuse = find(any(A(:,find(p.variabletype)),2));
            b(dontuse)=[];
            A(dontuse,:)=[];
            if size(A,1) > 0
                
                if p.options.bmibnb.verbose
                    disp('* -Trying to find bound on variables using SDP...');
                end
                
                A = A(:,v);
                n = length(v);
                gamma = sdpvar(1);
                sdpvar tau1 tau2
                if ~isinf(upper)
                    U = upper;
                else
                    U = sdpvar(1);
                end
                H1 = [gamma zeros(1,n);zeros(n,1) -eye(n)];
                H2 = [U - p.nonshiftedQP.f p.nonshiftedQP.c(v)'/2;p.nonshiftedQP.c(v)/2 -p.nonshiftedQP.Q(v,v)];
                H3 = [b -A]'*[b -A];
                if ~isinf(upper)
                    sol = optimize([[H1 - tau1*H2 - tau2*H3]>=0, tau1>=0],gamma,sdpsettings('verbose',p.options.verbose,'solver',p.solver.sdpsolver.tag));
                    if sol.problem == 0
                        r = sqrt(max(0,value(gamma)));
                        p.lb(v) = max(p.lb(v),-r);
                        p.ub(v) = min(p.ub(v),r);
                        disp('* - (managed to recover a radius on all variables)');
                    else
                        disp('* - (failed)');
                    end
                else
                    % OK, we don't have an upper bound. Save an optimizer to be
                    % used later
                    try
                        p.boundGenerator.optimizer = optimizer([[H1 - tau1*H2 - tau2*H3]>=0, tau1>=0],gamma,sdpsettings('solver',p.solver.sdpsolver.tag,'verbose',p.options.verbose),U,gamma);
                        disp('* - (not possible, will delay until upper bound available)');
                    catch
                        % No solver available?
                        p.boundGenerator.optimizer = [];
                    end
                end
            end
        end
    end
end