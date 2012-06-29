function varargout = gams2yalmip(fileName,filenameOut);
% GAMS2YALMIP Converts GAMS model to YALMIP
%
%   [F,h] = GAMS2YALMIP(gamsfile,yalmipfile) converts the GAMS model in
%   the file 'gamsfile' to a YALMIP mpodel
%
%   Input
%    GAMSFILE     : Char with filename for GAMS model
%    YALMIPFILE   : Char with filename for YALMIP model (optional)
%
%   Output
%    F            : LMI object with constraints  (optional)
%    h            : SDPVAR object with objective (optional)

% Author Based on original implementation by M. Kojima (readGMS.m)
% $Id: gams2yalmip.m,v 1.15 2007-08-06 07:50:21 joloef Exp $

writetofile = (nargout == 0) | (nargin>1);
fileName = [fileName '.gms'];
fileName = strrep(fileName,'.gms.gms','.gms');
% Reading the GAMs file "fileName"
fileIDX = fopen(fileName,'r');
if fileIDX==-1
     error('File not found')
end

if writetofile
    if nargin == 2
        filenameOut = [filenameOut '.m'];
        filenameOut = strrep(filenameOut,'.m.m','.m');
    else
        filenameOut = strrep(fileName,'.gms','.m');
    end
    fileOUT = fopen(filenameOut,'wb');
    if fileOUT==-1
        error('Could not create output-file')
    end
end

statusSW = 1;
noOfEquations = 0;
pp = 1;
posVarNames = [];
negVarNames = [];
lastline = '';
while statusSW == 1
    [statusSW,oneLine] = getOneLine(fileIDX);
    if statusSW == 1
        lastline = oneLine;
        [keyword,oneLine] = strtok(oneLine);
        if strcmp('Variables',keyword)
            varNames = [];
            p = 0;
            [varNames,p,moreSW] = getListOfNames(oneLine,varNames,p);
            while moreSW == 1
                [statusSW,oneLine] = getOneLine(fileIDX);
                [varNames,p,moreSW] = getListOfNames(oneLine,varNames,p);
            end
            noOfVariables = size(varNames,2);
            lbd = -inf* ones(1,noOfVariables);
            ubd = inf* ones(1,noOfVariables);
            fixed = lbd;
        elseif strcmp('Positive',keyword)
            [keyword,oneLine] = strtok(oneLine);
            if strcmp('Variables',keyword)
                p = 0;
                [posVarNames,p,moreSW] = getListOfNames(oneLine,posVarNames,p);
                while moreSW == 1
                    [statusSW,oneLine] = getOneLine(fileIDX);
                    [posVarNames,p,moreSW] = getListOfNames(oneLine,posVarNames,p);
                end
            end
        elseif strcmp('Negative',keyword)
            [keyword,oneLine] = strtok(oneLine);
            if strcmp('Variables',keyword)
                p = 0;
                [negVarNames,p,moreSW] = getListOfNames(oneLine,negVarNames,p);
                while moreSW == 1
                    [statusSW,oneLine] = getOneLine(fileIDX);
                    [negVarNames,p,moreSW] = getListOfNames(oneLine,negVarNames,p);
                end
            end
        elseif strcmp('Equations',keyword)
            equationNames  = [];
            p = 0;
            [equationNames,p,moreSW] = getListOfNames(oneLine,equationNames,p);
            while moreSW == 1
                [statusSW,oneLine] = getOneLine(fileIDX);
                [equationNames,p,moreSW] = getListOfNames(oneLine,equationNames,p);
            end
            noOfEquations = size(equationNames,2);
            listOfEquations = [];
        elseif pp <= noOfEquations
            if strcmp(strcat(equationNames{pp},'..'),keyword)
                oneLinetmp = oneLine;		% to remove blank around *
                while ~isempty(strfind(oneLinetmp,' *')) | ~isempty(strfind(oneLinetmp,'* '))
                    if ~isempty(strfind(oneLinetmp, ' *'))
                        loca = strfind(oneLinetmp,' *');
                        loca = loca -1;
                        oneLinetmp=strcat(oneLinetmp(1:loca),oneLinetmp(loca+2:size(oneLinetmp,2)));
                    elseif ~isempty(strfind(oneLinetmp, '* '))
                        loca = strfind(oneLinetmp,'* ');
                        oneLinetmp=strcat(oneLinetmp(1:loca),oneLinetmp(loca+2:size(oneLinetmp,2)));
                    end
                end
                oneLine = oneLinetmp;
                listOfEquations{pp} = oneLine;
                pp = pp+1;
            end
        elseif (0 < noOfEquations) & (noOfEquations < pp)
            goon = 1;
            while goon
                [oneVarName,bound] = strtok(keyword,'. ');
                for i=1:noOfVariables
                    if strcmp(oneVarName,varNames{i})
                        asciiVal = strtok(oneLine,' =;');
                        if strcmp(bound,'.lo') %| strcmp(bound,'.l')
                            lbd(1,i) = str2num(asciiVal);
                        elseif strcmp(bound,'.up') %| strcmp(bound,'.u')
                            ubd(1,i) = str2num(asciiVal);
                        elseif strcmp(bound,'.fx')
                            fixed(1,i) = str2num(asciiVal);
                        end
                    end
                end
                if strfind(oneLine,';')
                    oneLine = oneLine(min(strfind(oneLine,';'))+1:end);
                    [keyword,oneLine] = strtok(oneLine);
                    goon = ~isequal('',keyword);
                else
                    goon = 0;
                end
            end
        end
    end
end

% Figure out objective from last line
minimize = 1;
dirstart = strfind(lastline,'minimizing ');
obj = '[]';
if ~isempty(dirstart)
    [aux,obj] = strtok(lastline(dirstart:end));
else
    dirstart = strfind(lastline,'maximizing ');
    if ~isempty(dirstart)
        [aux,obj] = strtok(lastline(dirstart:end));
        minimize = -1;
    else
        minimize = 0;
    end
end
obj = strrep(strrep(obj,';',''),' ','');

objective_in_equations = 0;
for i = 1:length(listOfEquations)
    % If the objective variable is found in several equations, we define it
    % as a variable, and add all constraints, instead of assigninging it
    % from the typical objvar + f(x) =E= 0 expression
    if ~isempty(strfind(listOfEquations{i},obj))
        objective_in_equations = objective_in_equations  +1;
    end
end
if objective_in_equations>1
   treat_obj_as_var = 1;
else
   treat_obj_as_var = 0;  
end

if writetofile
    
    fprintf(fileOUT,['%% Model generated from ' fileName '\n']);
    [v,d] = yalmip('ver');
    fprintf(fileOUT,['%% Created ' datestr(now) ' using YALMIP R' d '\n\n']);
    fprintf(fileOUT,'%% Setup a clean YALMIP environment \n');
    fprintf(fileOUT,'yalmip(''clear'') \n\n');

    %fprintf(fileOUT,'%% Define non-standard operators \n');
    %fprintf(fileOUT,'sqr = @(x) x.*x;\n\n');

    % Define all variables, except objvar
    fprintf(fileOUT,'%% Define all variables \n');
end

for i = 1:length(varNames)
    eval([varNames{i} ' = sdpvar(1);']);
    if writetofile & (~isequal(varNames{i},obj) | treat_obj_as_var)
        fprintf(fileOUT,[varNames{i} ' = sdpvar(1);\n']);
    end
end
if writetofile
    fprintf(fileOUT,'\n');
end

if minimize
    if ~treat_obj_as_var
        if writetofile & objective_in_equations<=1
            fprintf(fileOUT,'%% Define objective function \n');
        end
        % find objvar + ... == 0
        for i = 1:length(listOfEquations)
            if strfind(listOfEquations{i},obj)
                %objeq = strrep(listOfEquations{i},'=E=','==');
                objeqL = listOfEquations{i}(1:strfind( listOfEquations{i},'=E=')-1);
                objeqR = listOfEquations{i}(strfind( listOfEquations{i},'=E=')+3:end);
                objeqR = strrep(objeqR,';','');

                % put objective on left side
                if strfind(objeqR,obj)
                    temp = objeqL;
                    objeqL = objeqR;
                    objeqR = objeqL;
                end
                k = strfind(objeqL,obj);
                prevplus = strfind(objeqL(1:k-1),'+');
                prevminus = strfind(objeqL(1:k-1),'-');
                if isempty(prevplus) & isempty(prevminus)
                    thesign = 1;
                else
                    prevsign = objeqL(max([prevplus prevminus]));
                    if isequal(prevsign,'+')
                        thesign = 1;
                    else
                        thesign = -1;
                    end
                end
                thesign = thesign*minimize;

                obj = [strrep(objeqL,obj,'0') '-( ' objeqR ')'];
                obj = strrep(obj,'**','^');
                obj = strrep(obj,'sqrt(','sqrtm(');
                obj = strrep(obj,'errorf(','erf(');
              %  obj = strrep(strrep(strrep(obj,'(  0)','0'),'0 -0','0'),'+ 0 ','');
              %  obj = strrep(strrep(strrep(obj,'(  0)','0'),'0 -0','0'),'+ 0)','');
                obj = strrep(obj,' + ','+');
                obj = strrep(obj,'+0)',')');
                obj = strrep(obj,'POWER','power');

                obj(obj==' ') = '';
                if writetofile
                    if thesign == -1
                        fprintf(fileOUT,['objective = ' obj ';\n']);
                    else
                        obj = ['-(' obj ')'];
                        obj = strrep(obj,'+0)',')');
                        obj = strrep(obj,'- ','-');
                        fprintf(fileOUT,['objective = ' obj ';\n']);
                    end
                end
                objsdp = (-thesign)*eval(obj);
                listOfEquations = {listOfEquations{1:i-1},listOfEquations{i+1:end}};
                break
            end
        end
        if writetofile
            fprintf(fileOUT,'\n');
        end
    end
end

% Convert to YALMIP syntax
for i = 1:length(listOfEquations)
    listOfEquations{i} = strrep(listOfEquations{i},'=E=','==');
    listOfEquations{i} = strrep(listOfEquations{i},'=L=','<=');
    listOfEquations{i} = strrep(listOfEquations{i},'=G=','>=');
end

% Add variable bounds
for i = 1:length(varNames)
    if any(strcmp(varNames{i},posVarNames))
        lbd(i) = 0;
    end
    if any(strcmp(varNames{i},negVarNames))
        ubd(i) = 0;
    end
    if ~isequal(varNames{i},obj) | treat_obj_as_var
        string = '';
        if ~isinf(fixed(i))
            string = [string num2str(fixed(i)) ' == ' varNames{i}];
        elseif ~isinf(lbd(i))
            string = [string num2str(lbd(i)) ' <= ' varNames{i}];
            if ~isinf(ubd(i))
                string = [string ' <= '  num2str(ubd(i))];
            end
        elseif ~isinf(ubd(i))
            string = [string varNames{i} ' <= ' num2str(ubd(i))  ];
        end
        if ~isequal(string,'')
            listOfEquations{end+1} = string;
        end
    end
end

if length(listOfEquations)>0
    if writetofile
        fprintf(fileOUT,'%% Define constraints \n');
    end
    F = set([]);
    if writetofile
        fprintf(fileOUT,['F = set([]);' '\n']);
    end
    for i = 1:length(listOfEquations)
        listOfEquations{i} = strrep(listOfEquations{i},';','');
        %        string = ['F = F + set(' listOfEquations{i} ',' '''' listOfEquations{i} '''' ');'];
        string = ['F = [F, ' strtrim(listOfEquations{i}) '];'];
        string = strrep(string,'**','^');
        string = strrep(string,'sqrt(','sqrtm(');
        string = strrep(string,'errorf(','erf(');
        string = strrep(string,'+',' + ');        
        string = strrep(string,'  ','');
        string = strrep(string,'  ','');
        string = strrep(string,' - ','-');        
        string = strrep(string,'POWER','power');            
        string = strtrim(string);
         string(string==' ') = '';
        eval(string);
        if writetofile
            fprintf(fileOUT,[string '\n']);
        end
    end
    if writetofile
        fprintf(fileOUT,'\n');
    end
else
    F = set([]);
    if writetofile
        fprintf(fileOUT,'%% Define constraints \n');
    end
    if writetofile
        fprintf(fileOUT,['F = set([]);' '\n']);
    end
end


if writetofile
    fprintf(fileOUT,'%% Solve problem\n');
    if treat_obj_as_var
        fprintf(fileOUT,'sol = solvesdp(F,objvar,sdpsettings(''solver'',''bmibnb'',''allownonconvex'',1));\n');
    else
        fprintf(fileOUT,'sol = solvesdp(F,objective,sdpsettings(''solver'',''bmibnb'',''allownonconvex'',1));\n');        
    end
    fprintf(fileOUT,'mbg_assertfalse(sol.problem)\n');        
    fprintf(fileOUT,'mbg_asserttolequal(double(objective), , 1e-2);');        

    fclose(fileOUT);
end

if nargout > 0
    varargout{1} = F;
    if nargout > 1
        varargout{2} = objsdp;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [statusSW,oneLine] = getOneLine(dataFID);
flowCTRL = 0;
oneLine = '';
while (feof(dataFID)==0) & (flowCTRL== 0)
    inputLine = fgetl(dataFID);
    %	inputLine
    len = length(inputLine);
    if (len > 0) & (inputLine(1)~='*')
        p=1;
        while (p<=len) & (inputLine(p)==' ')
            p = p+1;
        end
        if (p<=len) % & (inputLine(p) ~= '*')
            %			oneLine
            %			inputLine
            % Kojima 11/06/04; to meet MATLAB 5.2
            if isempty(oneLine)
                oneLine = inputLine(p:len);
            else
                oneLine = strcat(oneLine,inputLine(p:len));
            end
            % Kojima 11/06/04; to meet MATLAB 5.2
            %            temp = strfind(inputLine,';');
            temp = findstr(inputLine,';');
            if isempty(temp) == 0
                flowCTRL=1;
            end
        end
    end
end
if flowCTRL==0
    oneLine = '';
    statusSW = -1;
else
    statusSW = 1;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varNames,p,moreSW] = getListOfNames(oneLine,varNames,p);
while (length(oneLine) > 0)
    [oneName,remLine] = strtok(oneLine,' ,');
    if length(oneName) > 0
        p = p+1;
        varNames{p} = oneName;
    end
    oneLine = remLine;
end
lenLastVar = length(varNames{p});
if varNames{p}(lenLastVar) == ';'
    moreSW = 0;
    if lenLastVar == 1
        p = p-1;
    else
        varNames{p} = varNames{p}(1:lenLastVar-1);
    end
else
    moreSW = 1;
end
return
