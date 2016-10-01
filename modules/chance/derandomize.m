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
[randomVariables,map] = mergeDistributions(randomVariables);
groupedChanceConstraints = groupchanceconstraints(F);

[Fderandomized,eliminatedConstraints] = deriveChanceModel(groupedChanceConstraints,randomVariables);
Fderandomized = Fderandomized + F(find(keep)) + F(find(keep(~eliminatedConstraints)));
Fderandomized(find(is(Fderandomized,'random'))) = [];

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
                    case {'normal','normalm'}
                        theMean    = randomVariables{uncertaintyGroup}.distribution.parameters{2};
                        covariance = randomVariables{uncertaintyGroup}.distribution.parameters{3};
                        gamma = icdf('normal',confidencelevel,0,1);
                        if min(size(covariance))==1
                            covariance = diag(covariance);
                        end
                        if isa(covariance,'sdpvar')
                            error('Covariance cannot be an SDPVAR in normal distribution. Maybe you meant to use factorized covariance in ''normalf''');
                        end
                        if isempty(x)
                            Fderandomized = [Fderandomized, b + c_wTbase*theMean >= gamma*norm(chol(covariance)*(c_wTbase'))];
                        else
                            if isa(gamma,'sdpvar')
                                % This will be a nasty nonlinear model, so
                                % we cannot use SOCP-based norm operator
                                e = chol(covariance)*(AAA'*x+c_wTbase');
                                %Fderandomized = [Fderandomized, b + c_wTbase*theMean - (ccc + AAA*theMean)'*x >= gamma*sqrtm(e'*e)];
                                Fderandomized = [Fderandomized, b + c_wTbase*theMean - (ccc + AAA*theMean)'*x >= gamma*norm_callback(e)];
                            else
                                Fderandomized = [Fderandomized, b + c_wTbase*theMean - (ccc + AAA*theMean)'*x >= gamma*norm(chol(covariance)*(AAA'*x+c_wTbase'))];
                            end
                        end
                        eliminatedConstraints(ic)=1;
                    case 'normalf'
                        theMean    = randomVariables{uncertaintyGroup}.distribution.parameters{2};
                        covariance = randomVariables{uncertaintyGroup}.distribution.parameters{3};
                        gamma = icdf('normal',confidencelevel,0,1);
                        if isa(gamma,'sdpvar')
                            e = covariance*(AAA'*x+c_wTbase');
                            Fderandomized = [Fderandomized, b + c_wTbase*theMean - (ccc + AAA*theMean)'*x >= gamma*norm_callback(e)];
                        else
                            Fderandomized = [Fderandomized, b + c_wTbase*theMean - (ccc + AAA*theMean)'*x >= gamma*norm(covariance*(AAA'*x+c_wTbase'))];
                        end
                        eliminatedConstraints(ic)=1;
                    otherwise
                        error('Distribution not supported');
                end
            end
        end
    end
end

function [merged,map] = mergeDistributions(randomVariables);

merged = {};
map = [];
used = zeros(1,length(randomVariables));
mergedindex = 1;
for i = 1:length(randomVariables)
    if ~used(i)
        this = randomVariables{i};
        used(i) = 1;
        for j = i+1:length(randomVariables)
            if ~used(j)
                that = randomVariables{j};
                if strcmp(func2str(this.distribution.name),'random') && strcmp(func2str(that.distribution.name),'random')
                    if mergable(this.distribution.parameters{1},that.distribution.parameters{1})
                        % Same distribution                                                  
                        this = merge2distributions(this,that);                            
                        map(j) = mergedindex;
                        used(j) = 1;                      
                    end
                end
            end
        end
        merged{mergedindex} = this;
        map(i) = mergedindex;
        mergedindex = mergedindex + 1;
    end
end
    
    
function C = merge2distributions(A,B)

C = A;
C.variables = [A.variables;B.variables];
if any(strcmp(A.distribution.parameters{1},{'normal','normalf','normalm'})) || any(strcmp(A.distribution.parameters{1},{'normal','normalf','normalm'}))

    % A bit messy as there are two additional forms for normal
    % distributions, used to define mv normals, and normals with a
    % factorized covariance
    
    % First, normalize to matrix format
    dimAvar = size(A.distribution.parameters{3});
    dimBvar = size(B.distribution.parameters{3});
    if dimAvar(1) ~= dimAvar(2)
        varA = diag(A.distribution.parameters{3});
    else
        varA = A.distribution.parameters{3};
    end
    if dimBvar(1) ~= dimBvar(2)
        varB = diag(B.distribution.parameters{3});
    else
        varB = B.distribution.parameters{3};
    end
    % Is A in factor form, but not B, and vice versa. If so, put the other
    % one in factor form too    
    if strcmp(A.distribution.parameters{1},'normalf') &&  any(strcmp(B.distribution.parameters{1},{'normal','normalm'}))
        varB = chol(varB);        
    elseif strcmp(B.distribution.parameters{1},'normalf') &&  any(strcmp(A.distribution.parameters{1},{'normal','normalm'}))
        varA = chol(varA);        
    end
    varC = blkdiag(varA,varB);
    
    % Sort out various combinations of normal stuff in different forms
    if strcmp(A.distribution.parameters{1},'normal') && strcmp(B.distribution.parameters{1},'normal')
        varC = diag(varC);
        C.distribution.parameters{1} = 'normal';
        C.distribution.parameters{2} = [A.distribution.parameters{2};B.distribution.parameters{2}];
        C.distribution.parameters{3} = diag(varC);        
    elseif strcmp(A.distribution.parameters{1},'normalf') ||  any(strcmp(B.distribution.parameters{1},'normalf'))
        C.distribution.parameters{1} = 'normalf';
        C.distribution.parameters{2} = [A.distribution.parameters{2};B.distribution.parameters{2}];
        C.distribution.parameters{3} = varC;        
    else
        C.distribution.parameters{1} = 'normalm';
        C.distribution.parameters{2} = [A.distribution.parameters{2};B.distribution.parameters{2}];
        C.distribution.parameters{3} = varC;        
    end        
else
    for k = 2:length(A.distribution.parameters)
        C.distribution.parameters{k} = [A.distribution.parameters{k};B.distribution.parameters{k}];
    end    
end
    

function yesno = mergable(a,b)
if isequal(a,b)
    yesno = 1
elseif any(strcmp(a,{'normal','normalf','normalm'})) && any(strcmp(b,{'normal','normalf','normalm'}))
    yesno = 1;
else
    yesno = 0;
end