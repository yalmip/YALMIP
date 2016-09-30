function Fderandomized = derandomize(F)

% Goes through all probability constraints and checks for cases where we
% can use analytic expressions.

% Find chance constraints
chanceDeclarations = find(is(F,'chance'));
% Find variables with attached distributions
randomDeclarations = find(is(F,'random'));

if isempty(randomDeclarations)
    Fderandomized = F;
    return
end

keep = ones(length(F),1);
keep(randomDeclarations)=0;
keep(chanceDeclarations)=0;
randomVariables = extractRandomDefinitions(F(randomDeclarations));
groupedChanceConstraints = groupchanceconstraints(F);

[Fderandomized,eliminatedConstraints] = deriveChanceModel(groupedChanceConstraints,randomVariables);

Fderandomized = Fderandomized + F(find(keep)) + F(find(keep(~eliminatedConstraints))) + F(randomDeclarations);


function [Fderandomized,eliminatedConstraints] = deriveChanceModel(groupedChanceConstraints,randomVariables);

Fderandomized = [];
eliminatedConstraints = zeros(length(groupedChanceConstraints),1);

for uncertaintyGroup = 1:length(randomVariables)
    
    wVars = getvariables(randomVariables{uncertaintyGroup}.variables);
    
    for ic = 1:length(groupedChanceConstraints)
        if length(groupedChanceConstraints{ic})>1
            error('Joint chance constraint not supported');
        end
        if ~is(groupedChanceConstraints{ic},'elementwise')
            error('Only elementwise chance constraints supported')
        end
        X = sdpvar(groupedChanceConstraints{ic});
        if length(X)>1
            error('Only single elementwise chance constraints supported')
        end
        
        % OK, simple linear inequality
        allVars = depends(X);
        if ~isempty(intersect(wVars,allVars))
            xVars = setdiff(allVars,wVars);
            x = recover(xVars);
            w = recover(wVars);
            
            b = [];
            A = [];
            % Some pre-calc
            xw = [x;w];
            xind = find(ismembc(getvariables(xw),getvariables(x)));
            wind = find(ismembc(getvariables(xw),getvariables(w)));
            [Qs,cs,fs,dummy,nonquadratic] = vecquaddecomp(X,xw);
            c_wTbase = [];
            AAA = [];
            ccc = [];
            for i = 1:length(X)
                Q = Qs{i};
                c = cs{i};
                f = fs{i};
                if nonquadratic
                    error('Constraints can be at most quadratic, with the linear term uncertain');
                end
                Q_ww = Q(wind,wind);
                Q_xw = Q(xind,wind);
                Q_xx = Q(xind,xind);
                c_x = c(xind);
                c_w = c(wind);
                
                %b = [b;f + c_w'*w];
                %A = [A;-c_x'-w'*2*Q_xw'];
                % A = [A -c_x-2*Q_xw*w];
                AAA = [AAA;sparse(-2*Q_xw)];
                ccc = [ccc;-sparse(c_x)];
                b = [b;f];
                c_wTbase = [c_wTbase;c_w'];
            end
            % b = b + c_wTbase*w;
            % A = reshape(ccc + AAA*w,size(c_x,1),[]);
                                                                      
            if strcmp(func2str(randomVariables{uncertaintyGroup}.distribution.name),'random')                
                confidencelevel = struct(groupedChanceConstraints{ic}).clauses{1}.confidencelevel;
                switch randomVariables{uncertaintyGroup}.distribution.parameters{1}
                    case 'normal'
                        theMean    = randomVariables{uncertaintyGroup}.distribution.parameters{2};
                        covariance = randomVariables{uncertaintyGroup}.distribution.parameters{3};
                        gamma = icdf('normal',confidencelevel,0,1);
                        if isa(covariance,'sdpvar')
                            error('Covariance cannot be an SDPVAR in normal distribution. Maybe you meant to use factorized covariance in ''normalf''');
                        end
                        if isempty(x)
                            Fderandomized = [Fderandomized, b + c_wTbase*theMean >= gamma*norm(chol(covariance)*(c_wTbase'))];
                        else
                            Fderandomized = [Fderandomized, b + c_wTbase*theMean - (ccc + AAA*theMean)'*x >= gamma*norm(chol(covariance)*(AAA'*x+c_wTbase'))];
                        end
                        eliminatedConstraints(ic)=1;
                    case 'normalf'
                        theMean    = randomVariables{uncertaintyGroup}.distribution.parameters{2};
                        covariance = randomVariables{uncertaintyGroup}.distribution.parameters{3};
                        gamma = icdf('normal',confidencelevel,0,1);
                        Fderandomized = [Fderandomized, b + c_wTbase*theMean - (ccc + AAA*theMean)'*x >= gamma*norm(covariance*(AAA'*x+c_wTbase'))];
                        eliminatedConstraints(ic)=1;
                    otherwise
                        error('Distribution not supported');
                end
            end
        end
    end
end