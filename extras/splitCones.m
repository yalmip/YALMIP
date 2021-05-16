function [model,failure] = splitCones(model)

[model,failure] = splitSOCPCones(model);
if ~failure && any(model.K.s)
    [model,failure] = splitSDPCones(model);
end
if ~failure && any(model.K.e)
    [model,failure] = splitEXPCones(model);
end


function [model,failure] = splitSOCPCones(model)
% Find top of SOCP data
top = startofSOCPCone(model.K);
newF = [];
failure = 0;
if sum(model.K.q) > 0
    for i = 1:length(model.K.q)
        % c+A*expconestuff(x) + B*x in socp cone
        A = model.F_struc(top:top+model.K.q(i)-1,1+model.evalVariables);
        if nnz(A)>0
            if nnz(A(2:end-1,:))>0
                % Operator in norm expression not allowe
                failure = 1;
            else
                % Detect cone([1/2+f/2;x;1/2-f/2]) which is
                % what yalmip writes f >= x'*x as
                a = model.F_struc(top,:);
                b = model.F_struc(top+model.K.q(i)-1,:);
                c = a+b;
                if c(1)==1 && ~any(c(2:end))
                    % This is left-hand of f >= x'*x
                    d = a-b;
                    % rewrite as write as f >= t, t >= x'*x
                    model.F_struc(top,:)=0;
                    model.F_struc(top+model.K.q(i)-1,:)=0;
                    model.F_struc(top,1)=1/2;
                    model.F_struc(top+model.K.q(i)-1,1)=1/2;
                    model.F_struc(top,end+1)=1/2;
                    model.F_struc(top+model.K.q(i)-1,end)=-1/2;
                    d(1,end+1)=-1;
                    newF = addRow(newF,d);
                else
                    failure = 1;
                end
            end
        end
    end
end
if ~failure
    model = appendLPBlock(model,newF);
end

function [model,failure] = splitSDPCones(model)
% Find top of SDP data
top = startofSDPCone(model.K);
newF=[];
failure = 0;
if sum(model.K.s) > 0
    for i = 1:length(model.K.s)
        % Complete SDP basis
        B=model.F_struc(top:top+model.K.s(i)^2-1,1+model.evalVariables);
        if nnz(B)>0
            involvedRows = find(any(B,2));
            diagonalIndex = find(speye(model.K.s(i)));
            if any(~ismember(involvedRows,diagonalIndex))
                % Operator in non-diagonal
                failure = 1;
            else
                [~,loc] = ismember(involvedRows,diagonalIndex);
                for j = 1:length(loc)
                    % Detected operator in a diagonal
                    % Extract diagonal term a(x) and create new a(x) >= t
                    row = top + diagonalIndex(loc(j))-1;
                    a = model.F_struc(row,:);
                    a(1,end+1)=-1;
                    newF = addRow(newF,a);
                    % Replace diagonal with t
                    model.F_struc(row,:)=0;
                    model.F_struc(row,end+1) = 1;
                end
            end
        end
    end
end
if ~failure
    model = appendLPBlock(model,newF);
end

function [model,failure] = splitEXPCones(model)
% Find top of SDP data
top = startofEXPCone(model.K);
newF=[];
failure = 0;
if sum(model.K.e) > 0
    for i = 1:model.K.e       
        B=model.F_struc(top:top+2,1+model.evalVariables);
        % x3 >= x2 exp (x1/x2)
        % split x3 >= u, u >= x2 exp (v/x2), u >= x1
        if nnz(B)>0
            if any(B(1,:))
                a = -model.F_struc(top,:);
                a(1,end+1) = 1;
                newF = addRow(newF,a);
                % Replace v
                model.F_struc(top,:)=0;
                model.F_struc(top,end+1) = 1;                 
            end
            if any(B(3,:))
                a = model.F_struc(top+2,:);
                a(1,end+1) = -1;
                newF = addRow(newF,a);
                % Replace u
                model.F_struc(top+2,:)=0;
                model.F_struc(top+2,end+1) = 1;  
            end
            if any(B(2,:))
                failure = 1;
            end                                 
        end
        top = top + 3;
    end
end
if ~failure
    model = appendLPBlock(model,newF);
end

function Aa=addRow(A,a)
if isempty(A)
    Aa=a;
else
    A(end,length(a))=0;
    Aa=[A;a];
end

function model = appendLPBlock(model,newF)
model.F_struc = [model.F_struc(1:model.K.f+model.K.l,:);
    newF;
    model.F_struc(model.K.f+model.K.l+1:end,:)];
model.K.l = model.K.l + size(newF,1);
n = size(model.F_struc,2)-1;
model.c(n)=0;
model.Q(n,n)=0;