function diagnostic = callstrul(F,h,options)

F_new = [];

for i = 1:length(F)
    if ~is(F(i),'lmi')
        F_new = F_new + F(i);
    else
        X = sdpvar(F(i));
        [l,m,r]=factors(X);
        if isempty(m)
            F_new = F_new + F(i);
        else
            [L1,R1,A1,M1,negated_cont1,negated_disc1,epsilon1,delta1,numpos1,xindicies,Pindicies] = preprocess_constraint(X);
            F_new = F_new + assignschur(F(i),'HKM_schur_LR_structure',L1,R1,A1,M1,negated_cont1,negated_disc1,epsilon1,delta1,numpos1,xindicies,Pindicies);
        end
    end
end

if nargin < 2
    % Equalities are converted internally in SDPT3 to double-sided
    % inequalities. This messes up our Schur compiler, hence we do it
    % already outside SDPT3
    options = sdpsettings('solver','sdpt3','removeequalities',-1);
else
    options.solver = 'sdpt3'; 
    options.removeequalities = -1;
end
diagnostic = solvesdp(F_new,h,options);
