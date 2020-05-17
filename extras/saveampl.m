function solution = saveampl(varargin)
%SAVEAMPL Saves a problem definition in AMPL format
%
%    SAVEAMPL(F,h,'filename')    Saves the problem min(h(x)), F(x)>0 to the file filename
%    SAVEAMPL(F,h)               A "Save As"- box will be opened
%
% YALMIP is currently able to save problems with linear and non-linear
% element-wise inequality and equality constraints. Integer and binary
% variables are also supported.
%
% Note that YALMIP changes the variable names. Continuous variables
% are called x, binary are called y while z denotes integer variables.

F = varargin{1};
h = varargin{2};

% Expand nonlinear operators
options = sdpsettings;
[F2,failure,cause] = expandmodel(F,h,options);
if failure % Convexity propgation failed
    interfacedata = [];
    recoverdata = [];
    solver = '';
    diagnostic.solvertime = 0;
    diagnostic.problem = 14;
    diagnostic.info = yalmiperror(14,cause);
    return
end

%% FIXME: SYNC with expandmodel etc. Same in compileinterfacedata
setupBounds(F,options,yalmip('extvariables'));
solver.constraint.equalities.polynomial=0;
solver.constraint.binary=1;
solver.constraint.integer=0;
[F] = modelComplementarityConstraints(F,solver,[]);
        
        
%%
nvars = yalmip('nvars');
vars = depends(F);
vars = unique([vars depends(h)]);

binvars = yalmip('binvariables');
integervars = yalmip('intvariables');

for i = 1:length(F)
    if is(F(i),'binary')
        binvars = [binvars depends(F(i))];
    elseif  is(F(i),'integer')
        integervars = [integervars depends(F(i))];
    end
end

binvars = intersect(binvars,vars);
integervars = intersect(integervars,vars);

vars = setdiff(vars,union(integervars,binvars));
integervars = setdiff(integervars,binvars);
obj = amplexpr(h,vars,binvars,integervars);
constraints = {};

if ~isempty(F)
    for i = 1:length(F)
        if is(F(i),'element-wise')
            C = sdpvar(F(i));C=C(:);
            dummy = amplexpr(C,vars,binvars,integervars);
            for j = 1:length(C)
                if ~isempty(dummy{j})
                    constraints{end+1} = ['0 <= ' dummy{j}];
                end
            end           
        elseif is(F(i),'socp')
            C = sdpvar(F(i));C=C(:);
            dummy = amplexpr(C(1)^2-C(2:end)'*C(2:end),vars,binvars,integervars);
            constraints{end+1} = ['0 <= ' dummy{1}];
            dummy = amplexpr(C(1),vars,binvars,integervars);
            constraints{end+1} = ['0 <= ' dummy{1}];
            
        elseif is(F(i),'equality')
            C = sdpvar(F(i));C=C(:);
            dummy = amplexpr(C,vars,binvars,integervars);
            for j = 1:length(C)
                constraints{end+1} = ['0 == ' dummy{j}];
            end
        end
    end   
end

% Is a filename supplied
if nargin<3
    [filename, pathname] = uiputfile('*.mod', 'Save AMPL format file');
    if isa(filename,'double')
        return % User cancelled
    else
        % Did the user change the extension
        if isempty(strfind(filename,'.'))
            filename = [pathname filename '.mod'];
        else
            filename = [pathname filename];
        end
    end
else
    filename = varargin{3};
end

fid = fopen(filename,'w');
try

    %  fprintf(fid,['option randseed 0;\r\n']);
    if length(vars)>0
        fprintf(fid,['var x {1..%i};\r\n'],length(vars));
    end
    if length(binvars)>0
        fprintf(fid,['var y {1..%i} binary ;\r\n'],length(binvars));
    end
    if length(integervars)>0
        fprintf(fid,['var z {1..%i} integer ;\r\n'],length(integervars));
    end


    fprintf(fid,['minimize obj:  ' obj{1} ';'],max(vars));
    fprintf(fid,'\r\n');

    if length(constraints)>0
        for i = 1:length(constraints)
            constraints{i} = strrep(constraints{i},'mpower_internal','');
            fprintf(fid,['subject to constr%i: ' constraints{i} ';'],i);
            fprintf(fid,'\r\n');
        end
    end

    fprintf(fid,'solve;\r\n');
    if length(vars)>0
        fprintf(fid,'display x;\r\n');
    end
    if length(binvars)>0
        fprintf(fid,'display y;\r\n');
    end
    if length(integervars)>0
        fprintf(fid,'display z;\r\n');
    end

    fprintf(fid,'display obj;\r\n');
catch
    fclose(fid);
end
fclose(fid);


function symb_pvec = amplexpr(pvec,vars,binvars,integervars)

extVariables = yalmip('extvariables');
for pi = 1:size(pvec,1)
    for pj = 1:size(pvec,2)
        p = pvec(pi,pj);

        if isa(p,'double')
            symb_p = num2str(p,12);
        elseif isinf(getbasematrix(p,0))
            symb_p = [];
        else
            LinearVariables = depends(p);
            x = recover(LinearVariables);
            exponent_p = full(exponents(p,x));
            names = cell(length(LinearVariables),1);
            for i = 1:length(LinearVariables)               
                v1 = find(vars==LinearVariables(i));
                if ~isempty(v1)
                    names{i}=['x[' num2str(find(vars==LinearVariables(i))) ']'];
                else
                    v1 = find(binvars==LinearVariables(i));
                    if ~isempty(v1)
                        names{i}=['y[' num2str(find(binvars==LinearVariables(i))) ']'];
                    else
                        names{i}=['z[' num2str(find(integervars==LinearVariables(i))) ']'];
                    end
                end
            end
            for i = 1:length(LinearVariables)                
                v1 = find(extVariables==LinearVariables(i));
                if ~isempty(v1)
                    e = yalmip('extstruct',extVariables(v1));
                    inner = amplexpr(e.arg{1},vars,binvars,integervars);
                    names{i} = [e.fcn '(' inner{1} ')'];   
                    names{i} = strrep(names{i},'mpower_internal','');
                end
            end

            symb_p = '';
            if all(exponent_p(1,:)==0)
                symb_p = num2str(full(getbasematrix(p,0)),12);
                exponent_p = exponent_p(2:end,:);
            end

            for i = 1:size(exponent_p,1)
                coeff = getbasematrixwithoutcheck(p,i);
                switch full(coeff)
                    case 1
                        coeff='+';
                    case -1
                        coeff = '-';
                    otherwise
                        if coeff >0
                            coeff = ['+' num2str2(coeff)];
                        else
                            coeff=[num2str2(coeff)];
                        end
                end
                if strcmp(symb_p,'') & (strcmp(coeff,'+') | strcmp(coeff,'-'))
                    symb_p = [symb_p coeff symbmonom(names,exponent_p(i,:))];
                else
                    symb_p = [symb_p coeff '*' symbmonom(names,exponent_p(i,:))];
                end
            end
            if symb_p(1)=='+'
                symb_p = symb_p(2:end);
            end
        end
    if ~isempty(symb_p)
        symb_p = strrep(symb_p,'+*','+');
        symb_p = strrep(symb_p,'-*','-');
    end
        symb_pvec{pi,pj} = symb_p;
    end
end

function s = symbmonom(names,monom)
s = '';
for j = 1:length(monom)
    if monom(j)
        if strcmp(s,'')
            s = [s names{j}];
        else
            s = [s '*' names{j}];
        end
        
        if monom(j)~=1
            s = [s '^' num2str(monom(j))];
        end
    end
end

function s = num2str2(x)
s = num2str(full(x),12);
if isequal(s,'1')
    s = '';
end
if isequal(s,'-1')
    s = '-';
end
  
