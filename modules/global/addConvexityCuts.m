function pcut = addConvexityCuts(p)

% Quick hack for a particular example...
pcut = p;
if ~isempty(p.evalMap)
    pcut = emptyNumericalModel;
    for i = find(p.variabletype)
        involved = find(p.monomtable(i,:));
        if length(involved) == 2
            evals = find(ismember(involved,p.evalVariables));
            if length(evals) == 1
                theOperator =  involved(evals);
                theOperatorIndex = find(theOperator == p.evalVariables);
                switch p.evalMap{theOperatorIndex}.properties.convexity
                    case 'convex'
                        theOther = setdiff(involved,theOperator);
                        if p.variabletype(theOther) == 0 & isequal(p.evalMap{1}.variableIndex,theOther)
                            % x*f(x)
                            if strcmp(p.evalMap{theOperatorIndex}.properties.monotonicity,'increasing')
                                xM = (p.lb(theOther) + p.ub(theOther))/2;
                                fM = xM*feval(p.evalMap{theOperatorIndex}.fcn,xM);
                                dfM = p.evalMap{theOperatorIndex}.properties.derivative(xM) + xM*p.evalMap{theOperatorIndex}.properties.derivative(xM);
                                row = spalloc(1,1+length(p.variabletype),3);
                                % y >= f + df*(x-xm)
                                row(1) = -fM+dfM*xM;
                                row(1+i) = 1;
                                row(1+theOther) = -dfM;
                                pcut.F_struc = [ pcut.F_struc;row];
                                pcut.K.l = pcut.K.l + 1;
                            end
                        end
                    otherwise
                end
            end
        end
    end
    pcut = mergeNumericalModels(p,pcut);
end

