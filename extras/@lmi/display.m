function sys = display(X)
%display           Displays a SET object.

% Author Johan Löfberg
% $Id: display.m,v 1.12 2009-05-29 08:05:12 joloef Exp $

nlmi = size(X.clauses,2);
nlmi = length(X.LMIid);

if (nlmi == 0)
    disp('empty SET')
    return
end

lmiinfo{1} = 'Matrix inequality';
lmiinfo{2} = 'Element-wise inequality';
lmiinfo{3} = 'Equality constraint';
lmiinfo{4} = 'Second order cone constraint';
lmiinfo{5} = 'Rotated Lorentz constraint';
lmiinfo{7} = 'Integrality constraint';
lmiinfo{8} = 'Binary constraint';
lmiinfo{9} = 'KYP constraint';
lmiinfo{10}= 'Eigenvalue constraint';
lmiinfo{11}= 'Sum-of-square constraint';
lmiinfo{12}= 'Logic constraint';
lmiinfo{13}= 'Parametric declaration';
lmiinfo{14}= 'Low-rank data declaration';
lmiinfo{15}= 'Uncertainty declaration';
lmiinfo{16}= 'Distribution declaration';
lmiinfo{20}= 'Power cone constraint';
lmiinfo{30}= 'User defined compilation';
lmiinfo{40}= 'Generalized KYP constraint';
lmiinfo{50}= 'Special ordered set of type 2';
lmiinfo{51}= 'Special ordered set of type 1';
lmiinfo{52}= 'Semi-continuous variable';
lmiinfo{53}= 'Semi-integer variable';
lmiinfo{54} = 'Vectorized second order cone constraints';
lmiinfo{55}= 'Complementarity constraint';
lmiinfo{56}= 'Meta constraint';

headers = {'ID','Constraint','Type','Tag'};
rankVariables = yalmip('rankvariables');
extVariables = yalmip('extvariables');
if nlmi>0
    for i = 1:nlmi
        
        data{i,1} = ['#' num2str(i)];
        data{i,2} = X.clauses{i}.symbolic;
        data{i,3} = lmiinfo{X.clauses{i}.type};
        data{i,4} = '';
        if length(getvariables(X.clauses{i}.data)) == 1
            if any(ismember(getvariables(X.clauses{i}.data),rankVariables))
                 data{i,3} = 'Rank constraint';                
            end
        end
        
        if X.clauses{i}.type == 14
            
        elseif X.clauses{i}.type == 56
            data{i,3} = [data{i,3} ' (' X.clauses{i}.data{1} ')'];  
            data{i,4} = X.clauses{i}.handle;
              
        else
            classification = '';
           
            members = ismembcYALMIP(getvariables(X.clauses{i}.data),yalmip('intvariables'));
            if any(members)
                classification = [classification ',integer'];
            end

            if size(X.clauses{i},2)>1
                classification = [classification ',logic'];                
            end
           
            if ~isempty(X.clauses{i}.confidencelevel)            
                classification = [classification ',chance'];                
            end
           
            linearbilinearquadraticsigmonial = is(X.clauses{i}.data,'LBQS');
            if ~linearbilinearquadraticsigmonial(1)
                if linearbilinearquadraticsigmonial(4)
                    classification = [classification ',sigmonial'];
                elseif linearbilinearquadraticsigmonial(2)
                    classification = [classification ',bilinear'];
                elseif linearbilinearquadraticsigmonial(3)
                    classification = [classification ',quadratic'];
                else
                    classification = [classification ',polynomial'];
                end                
            end
            
            data{i,4} = X.clauses{i}.handle;
            if ~isreal(X.clauses{i}.data)                
                classification = [classification ',complex'];
            end
            %if ~isempty(intersect(getvariables(X.clauses{i}.data),extVariables))            
            members = ismembcYALMIP(getvariables(X.clauses{i}.data),extVariables);          
            if any(members)
                classification = [classification ',derived'];
            end
            
            if length(classification)==0
            else                
                data{i,3} = [data{i,3} ' (' classification(2:end) ')'];
            end

            if ismember(X.clauses{i}.type,[1 2 3 4 5 9]);
                data{i,3} = [data{i,3} ' ' num2str(size(X.clauses{i}.data,1)) 'x' num2str(size(X.clauses{i}.data,2))];
            end
        end

    end
end

% If no tags, don't show...
if length([data{:,4}])==0
    headers = {headers{:,1:3}};
    data = reshape({data{:,1:3}},length({data{:,1:3}})/3,3);
end



yalmiptable('',headers,data)

function x= truncstring(x,n)
if length(x) > n
    x = [x(1:n-3) '...'];
end

function x = fillstring(x,n)
x = [x blanks(n-length(x))];

