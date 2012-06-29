function [Nbegcol,Nlencol,Nrowndx,Nobjcnt,Nobjndx,cJacobian] = jacSparsity(interfacedata)

linear = setdiff(find(interfacedata.variabletype == 0),interfacedata.evalVariables);
oJacobian = zeros(length(linear),1);
variabletype = interfacedata.variabletype;
c = interfacedata.c;
F_struc = interfacedata.F_struc;
m = size(interfacedata.F_struc,1);

nonlinear = variabletype;
nonlinear(interfacedata.evalVariables) = 1;

for i = find(c(:)')%1:length(c)
%    if c(i)
        if nonlinear(i)%variabletype(i) | ismember(i,interfacedata.evalVariables)
            if ismember(i,interfacedata.evalVariables)
                j = find(i == interfacedata.evalVariables);
                variables = interfacedata.evalMap{j}.variableIndex;
            else
                variables = find(interfacedata.monomtable(i,:));
            end
            oJacobian(find(ismember(linear,variables))) = 1;
        end
%    end
end
if m > 0
    cJacobian = zeros(m,length(linear));
    for i = 1:size(F_struc,2)-1
        %if nonlinear(i)%variabletype(i) | ismembc(i,interfacedata.evalVariables)
            f = F_struc(:,i+1);
            for j = find(f(:)')%1:size(F_struc,1)
                %   if f(j)%F_struc(j,i+1)
                %                if variabletype(i) | ismembc(i,interfacedata.evalVariables)
                if ismembc(i,interfacedata.evalVariables)
                    k = find(i == interfacedata.evalVariables);
                    variables = interfacedata.evalMap{k}.variableIndex;
                elseif nonlinear(i)
                    variables = find(interfacedata.monomtable(i,:));
                else
                    variables = i;;
                end
                cJacobian(j,find(ismembc(linear,variables))) = 1;
            end            
       % else
            
       % end
    end 
    [Nbegcol,Nlencol,Nrowndx] = lindosparse(cJacobian);
else
    cJacobian = [];
    Nbegcol = [];
    Nrowndx = [];
    Nlencol = [];
    Nbegcol = [Nbegcol sum(Nlencol)];
end
oJacobian = double(oJacobian | any(interfacedata.Q(linear,linear),2));

Nobjndx = find(oJacobian) - 1;
Nobjcnt = length(Nobjndx);
if  isempty(Nobjndx)
    Nobjndx = [];
end


