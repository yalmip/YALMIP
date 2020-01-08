function filename = savecplexlp(varargin)
%SAVCPLEXLP Saves a problem definition in CPLEX-LP format
%
%    SAVCPLEXLP(F,h,'filename')    Saves the problem min(h(x)), F(x)>0 to the file filename
%    SAVCPLEXLP(F,h)               A "Save As"- box will be opened
%

F = varargin{1};
h = varargin{2};

[aux1,aux2,aux3,model] = export(F,h);

% Check so that it really is an LP
% if any(any(model.Q)) | any(model.variabletype)
% Only check the variabletype
if any(model.variabletype)
    error('This is not an LP or QP');
end

c = model.c;
Q = 2*full(model.Q);
b = model.F_struc(:,1);
A = -model.F_struc(:,2:end);
if model.K.f>0
    Aeq = A(1:model.K.f,:);
    beq = b(1:model.K.f,:);
    A(1:model.K.f,:) = [];
    b(1:model.K.f) = [];
else
    Aeq = [];
    beq = [];
end

lb = model.lb;
ub = model.ub;

[lb,ub,A,b] = remove_bounds_from_Ab(A,b,lb,ub);
[lb,ub,Aeq,beq] = remove_bounds_from_Aeqbeq(Aeq,beq,lb,ub);

% Is a filename supplied
if nargin<3
    [filename, pathname] = uiputfile('*.lp', 'Save LP format file');
    if isa(filename,'double')
        return % User cancelled
    else
        % Did the user change the extension
        if isempty(findstr(filename,'.'))
            filename = [pathname filename '.lp'];
        else
            filename = [pathname filename];
        end
    end
else
    filename = varargin{3};
end

fid = fopen(filename,'w');

obj = strrep(lptext(c(:)'),'+ -','-');
    
fprintf(fid,['Minimize\r\n obj:  ' obj(1:end-2) '']);
if any(any(Q))
    obj = qptext(Q);
    fprintf(fid,[' + [' obj(1:end-2) '] / 2']);
end
fprintf(fid,'\r\n');

fprintf(fid,['\r\n']);
fprintf(fid,['Subject To\r\n']);

for i = 1:length(b)
    rowtext = lptext(-A(i,:));
    rhs = sprintf('%0.20g',full(-b(i)));
    rowtext = [rowtext(1:end-2) ' >= ' rhs]; 
    fprintf(fid,[' c%i: ' strrep(rowtext,'+ -','-') ''],i);
    fprintf(fid,'\r\n');
end
for i = 1:length(beq)
    rowtext = lptext(-Aeq(i,:));
    rowtext = [rowtext(1:end-2) '== ' sprintf('%0.20g',full(-beq(i)))];    
    fprintf(fid,[' eq%i: ' strrep(rowtext,'+ -','-') ''],i);
    fprintf(fid,'\r\n');
end

if length(c)>length(model.binary_variables)
    fprintf(fid,['\r\nBounds\r\n']);
    for i = 1:length(c)
        %        if ~ismember(i,model.binary_variables)
        if isinf(lb(i)) & isinf(ub(i))
            fprintf(fid,[' x%i free\n\r'],i);
        elseif lb(i)==0 & isinf(ub(i))
            % Standard non-negative variable
        elseif isinf(ub(i))
            s = strrep(sprintf(['%0.20g <= x%i \r\n'],[lb(i) i ]),'Inf','inf');
            fprintf(fid,s);
        else
            s = strrep(sprintf(['%0.20g <= x%i <= %0.20g \r\n'],[lb(i) i ub(i)]),'Inf','inf');
            fprintf(fid,s);
        end
        %        end
    end
end

if length(model.binary_variables)>0
    fprintf(fid,['\r\n']);
    fprintf(fid,['Binary\r\n']);
    for i = 1:length(model.binary_variables)
        fprintf(fid,[' x%i\r\n'],model.binary_variables(i));
    end
end

if length(model.integer_variables)>0
    fprintf(fid,['\r\n']);
    fprintf(fid,['Integer\r\n']);
    for i = 1:length(model.integer_variables)
        fprintf(fid,[' x%i\r\n'],model.integer_variables(i));
    end
end

fprintf(fid,['\r\nEnd']);
fclose(fid);

function rowtext = lptext(a)
[aux,poss,vals] = find(a);
rowtext = sprintf('%0.20g x%d + ',reshape([vals(:) poss(:)]',[],1));
%rowtext = strrep(rowtext,'+ -','- ');
%rowtext(isspace(rowtext))=[];
%rowtext = strrep(rowtext,'+-','-');
%rowtext = strrep(rowtext,'-1x','-x');
%rowtext = strrep(rowtext,'+1x','+x');

function rowtext = qptext(Q)
n=size(Q,2);
q = diag(Q);
if any(q)
    i = find(q);
    rowtext = sprintf('%0.20g x%d ^2 + ',reshape([q(i) i]',[],1));
end
for i=1:n
    for j=i+1:n
        if ~(Q(i,j)+Q(j,i)==0)
        rowtext = [rowtext sprintf('%0.20g x%d * x%d + ',Q(i,j)+Q(j,i),i,j)];
        end
    end
end
rowtext = strrep(rowtext,'+ -','- ');