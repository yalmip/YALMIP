function display(X)
%DISPLAY (overloaded)

% Author Johan Löfberg
% $Id: display.m,v 1.2 2006-08-10 23:02:56 joloef Exp $

switch(X.typeflag)
    case {0,9}
        n = X.dim(1);
        m = X.dim(2);
        if (n*m==1)

            linearbilinearquadraticsigmonial = is(X,'LBQS');            
            if linearbilinearquadraticsigmonial(1)
                classification = 'Noncommuting linear scalar ';
            elseif linearbilinearquadraticsigmonial(4)
                classification = 'Noncommuting sigmonial scalar ';
            elseif linearbilinearquadraticsigmonial(2)
                classification = 'Noncommuting bilinear scalar ';
            elseif linearbilinearquadraticsigmonial(3)
                classification = 'Noncommuting quadratic scalar ';
            else
                classification = 'Noncommuting polynomial scalar ';
            end

            if ~isreal(X.basis)
                classification = [classification '(complex'];
            else
                classification = [classification '(real'];
            end

            if is(X,'compound')
                classification = [classification ', derived'];
            end

            if ~islinear(X)
                variables = getvariables(X);
                monomtable = yalmip('monomtable');
                if ((nnz(getbasematrix(X,0))==0) ) & (sum(diff(sum(monomtable(variables,:),2)))==0)
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
            if ~isnan(double(X))
                classification = [classification ', current value : ' num2str(double(X))];
            end
            classification = [classification ')'];

            disp([classification]);
        else

            if islinear(X)
                classification = 'Noncommuting linear matrix variable ';
            elseif is(X,'sigmonial')
                classification = 'Noncommuting sigmonial matrix variable ';
            elseif is(X,'bilinear')
                classification = 'Noncommuting bilinear matrix variable ';
            elseif is(X,'quadratic')
                classification = 'Noncommuting quadratic matrix variable ';
            else
                classification = 'Noncommuting polynomial matrix variable ';
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

            if is(X,'compound')
                info = [info ', derived'];
            end

            if X.typeflag==9
                info = [info ', KYP'];
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
            
            if (n<100) & (n==m)
                x = recover(xvars);
                if ~any(any(isnan(double(x))))
                    doubleX = double(X);
                    try
                    eigX = eig(doubleX);
                    info = [info ', eigenvalues between [' num2str(min(eigX)) ',' num2str(max(eigX)) ']'];               
                    catch
                    end
                end
            elseif n~=m
                x = recover(xvars);
                if ~any(any(isnan(double(x))))
                    doubleX = double(X);
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
    otherwise
end