function display(X)
%DISPLAY (overloaded)

switch(X.typeflag)
    case {0,9,40}
        n = X.dim(1);
        m = X.dim(2);
        if (n*m==1)

            vars = depends(X);
            if any(ismember(vars,yalmip('extvariables')))
                classification = 'Nonlinear scalar ';
            else                
                linearbilinearquadraticsigmonial = is(X,'LBQS');
                if linearbilinearquadraticsigmonial(1)
                    classification = 'Linear scalar ';
                elseif linearbilinearquadraticsigmonial(4)
                    classification = 'Signomial scalar ';
                elseif linearbilinearquadraticsigmonial(2)
                    classification = 'Bilinear scalar ';
                elseif linearbilinearquadraticsigmonial(3)
                    classification = 'Quadratic scalar ';
                else
                    classification = 'Polynomial scalar ';
                end
            end

            if ~isreal(X.basis)
                classification = [classification '(complex'];
            else
                classification = [classification '(real'];
            end

            if is(X,'compound')               
                if ~isequal(X.extra.opname,'')
                    classification = [classification ', models ''' X.extra.opname ''''];
                end
            end

            if ~islinear(X)
                variables = getvariables(X);
                monomtable = yalmip('monomtable');
                s = sum((monomtable(variables,:)),2);
               % s = sum(full(monomtable(variables,:)),2);
                if ((nnz(getbasematrix(X,0))==0) ) && all(s==s(1))
                    classification = [classification ', homogeneous'];
                end
            end
            if any(ismember(depends(X),yalmip('intvariables')))
                classification = [classification ', integer'];
            else
                if any(ismember(depends(X),yalmip('binvariables')))
                    classification = [classification ', binary'];
                else
                    if any(ismember(depends(X),yalmip('uncvariables')))
                        classification = [classification ', uncertain'];
                    else
                        if any(ismember(depends(X),yalmip('semicontvariables')))
                        classification = [classification ', semi-continuous'];   
                        end
                    end
                end
            end

            nvars = length(depends(X));
            if nvars == 1
                classification = [classification ', ' num2str(nvars) ' variable'];
            else
                classification = [classification ', ' num2str(nvars) ' variables'];
            end

            if any(ismember(depends(X),yalmip('parvariables')))
                classification = [classification ', parametric'];
            end
            if ~isnan(value(X))
                classification = [classification ', current value : ' num2str(value(X))];
            end
            classification = [classification ')'];

            disp([classification]);
        else

            vars = depends(X);
            if any(ismember(vars,yalmip('extvariables')))
                classification = 'Nonlinear matrix variable ';
            else
                if islinear(X)
                    classification = 'Linear matrix variable ';
                elseif is(X,'sigmonial')
                    classification = 'Signomial matrix variable ';
                elseif is(X,'bilinear')
                    classification = 'Bilinear matrix variable ';
                elseif is(X,'quadratic')
                    classification = 'Quadratic matrix variable ';
                else
                    classification = 'Polynomial matrix variable ';
                end
            end

            if isreal(X.basis)
                if issymmetric(X)
                    info = ' (symmetric, ';
                else
                    info = ' (full, ';
                end;
                info = [info 'real'];
            else
                if issymmetric(X)
                    info = ' (symmetric, ';
                elseif ishermitian(X)
                    info = ' (hermitian, ';
                else
                    info = ' (full, ';
                end;
                info = [info 'complex'];
            end;

            if X.typeflag==9
                info = [info ', KYP'];
            end
            
            if X.typeflag==40
                info = [info ', Generalized KYP'];
            end            
            
            if any(ismember(depends(X),yalmip('intvariables')))
                info = [info ', integer'];
            else
                if any(ismember(depends(X),yalmip('binvariables')))
                    info = [info ', binary'];
                else
                    if any(ismember(depends(X),yalmip('uncvariables')))
                        info = [info ', uncertain'];
                    end
                end
            end
            
            if any(ismember(depends(X),yalmip('parvariables')))
                info = [info ', parametric'];
            end
            
            xvars = depends(X);
            nvars = length(xvars);
            if nvars == 1
                info = [info ', ' num2str(nvars) ' variable'];
            else
                info = [info ', ' num2str(nvars) ' variables'];
            end
            
            n = X.dim(1);
            m = X.dim(2);
            
            if (n<100) && (n==m)
                x = recover(xvars);
                if ~any(any(isnan(value(x))))
                    doubleX = double(X);
                    try
                    eigX = eig(doubleX);
                    info = [info ', eigenvalues between [' num2str(min(eigX)) ',' num2str(max(eigX)) ']'];               
                    catch
                    end
                end
            elseif n~=m
                x = recover(xvars);
                if ~any(any(isnan(value(x))))
                    doubleX = value(X);
                    try                       
                        info = [info ', values in range [' num2str(min(min(doubleX))) ',' num2str(max(max(doubleX))) ']'];
                    catch
                    end
                end
            end
          
            info = [info ')'];
            disp([classification num2str(n) 'x' num2str(m) info]);
        end;
    case 1
        disp('Relational object');
    case 3
        disp('Relational object');
    case 4
        disp('Relational object');
    case 5
        disp('Cone object');
    case 6
        disp('logdet object');
    case 7
        disp('Integrality constraint');
    case 10
        disp('Eigenvalue object');
    case 11
        disp('SOS object');
    case 20
        disp('Power cone object');
        
    otherwise
end