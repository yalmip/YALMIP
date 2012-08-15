function Fderandomized = derandomize(F)

chanceDeclarations = find(is(F,'chance'));
randomDeclarations = find(is(F,'random'));
if isempty(randomDeclarations)
    error('Cannot derandomize, cannot find declaration of random variables');
end

keep = ones(length(F),1);
keep(randomDeclarations)=0;
keep(chanceDeclarations)=0;
randomVariables = extractRandomDefinitions(F(randomDeclarations));
groupedChanceConstraints = groupchanceconstraints(F);

Fderandomized = deriveChanceModel(groupedChanceConstraints,randomVariables)

Fderandomized = Fderandomized + F(find(keep));


function Fderandomized = deriveChanceModel(groupedChanceConstraints,randomVariables);

Fderandomized = [];
for ic = 1:length(groupedChanceConstraints)
    if length(groupedChanceConstraints{ic})>1
        error('Joint chance still not supported');
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
    wVars = [];for j = 1:length(randomVariables);wVars = [wVars getvariables(randomVariables{j}.variables)];end
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
 
   
   j = 1;
   confidencelevel = struct(groupedChanceConstraints{ic}).clauses{1}.confidencelevel;
   theMean    = randomVariables{j}.distribution.parameters{1};
   covariance = randomVariables{j}.distribution.parameters{2};
   gamma = icdf('normal',confidencelevel,0,1);
   
   switch randomVariables{j}.distribution.name
       case 'normal'
           theMean    = randomVariables{j}.distribution.parameters{1};
           covariance = randomVariables{j}.distribution.parameters{2};
           gamma = icdf('normal',confidencelevel,0,1);
           if isa(covariance,'sdpvar')
               error('Covariance cannot be an SDPVAR in normal distribution. Maybe you meant to use factorized covariance in ''normalf''');
           end
           Fderandomized = [Fderandomized, b + c_wTbase*theMean - (ccc + AAA*theMean)'*x >= gamma*norm(chol(covariance)*(AAA'*x+c_wTbase'))];
       case 'normalf'
           theMean    = randomVariables{j}.distribution.parameters{1};
           covariance = randomVariables{j}.distribution.parameters{2};
           gamma = icdf('normal',confidencelevel,0,1);           
           Fderandomized = [Fderandomized, b + c_wTbase*theMean - (ccc + AAA*theMean)'*x >= gamma*norm(covariance*(AAA'*x+c_wTbase'))];
       otherwise
           error('Distribution not supported');
   end
   
    
    
    
    
end