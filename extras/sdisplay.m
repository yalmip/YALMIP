function symb_pvec = sdisplay(pvec,precision)
%SDISPLAY Symbolic display of SDPVAR expression
%
% As variables in YALMIP do not have any actual names internally, the
% command is not guaranteed to recover the variable names you use in the
% work-space correctly. It tries to figure out the names of the variables
% by searching in the work-space for variables involving the same variable
% indicies, but this can fails in many instances.
%
% EXAMPLES
%  sdpvar x y
%  sdisplay(x^2+y^2)
%    ans =
%       'x^2+y^2'
%
% Optionally, we can supply an option to round the data to N decimal digits
%
%  sdisplay(pi*x1^2,4)
%    ans = '3.1416*x^2'

if nargin < 2
    precision = inf;
end

% Support displaying non-sdpvar input
if ~isa(pvec,'sdpvar')
    for r1=1:size(pvec,1)
        for r2=1:size(pvec,2)
            p = pvec(r1,r2);
            if isa(p,'double')
                symb_pvec{r1,r2} = num2str2(p,precision);
            else
                symb_pvec{r1,r2} = p;
            end
        end
    end
    if prod(size(symb_pvec))==1 & nargout==0
        display(symb_pvec{1,1});
        clear symb_pvec
    end
    return;
end

% First, some boooring stuff. we need to
% figure out the symbolic names and connect
% these names to YALMIPs variable indicies
W = evalin('caller','whos');

% First, sort the variables available in the work-space based
% on the creation. Early creation means it is more likely the
% relevant variable to display
createTime = [];
for i = 1:size(W,1)
    if strcmp(W(i).class,'sdpvar') || strcmp(W(i).class,'ncvar')
        keep(i) = 1;
        z = evalin('caller',['struct(' W(i).name ').extra.createTime;']);
        createTime = [createTime z];
    else
        keep(i) = 0;
    end
end
W = W(find(keep));
[sorted,index] = sort(createTime);
W = W(index);

global_LinearVariables = depends(pvec);
global_x = recover(global_LinearVariables);
[global_exponent_p,global_ordered_list] = exponents(pvec,global_x);
global_exponent_p = full(global_exponent_p);
global_names = cell(length(global_x),1);

for i = 1:size(W,1)
    if strcmp(W(i).class,'sdpvar') || strcmp(W(i).class,'ncvar')
        % Get the SDPVAR variable
        thevars = evalin('caller',W(i).name);
        
        % Distinguish 4 cases
        % 1: Sclalar varible x
        % 2: Vector variable x(i)
        % 3: Matrix variable x(i,j)
        % 4: Variable not really defined
        if is(thevars,'scalar') && is(thevars,'linear') && length(getvariables(thevars))==1 & isequal(getbase(thevars),[0 1])
            index_in_p = find(ismember(global_LinearVariables,getvariables(thevars)));
            
            if ~isempty(index_in_p)
                already = ~isempty(global_names{index_in_p});
                if already
                    already = (isempty(strfind(global_names{index_in_p},'internal')) | isempty(strfind(global_names{index_in_p},'ans')));
                end
            else
                already = 0;
            end
            
            if ~isempty(index_in_p) & ~already
                % Case 1
                global_names{index_in_p}=W(i).name;
            end
        elseif is(thevars,'lpcone')
            
            if size(thevars,1)==size(thevars,2)
                % Case 2
                vars = getvariables(thevars);
                indicies = find(ismember(vars,global_LinearVariables));
                for ii = indicies
                    index_in_p = find(ismember(global_LinearVariables,vars(ii)));
                    
                    if ~isempty(index_in_p)
                        already = ~isempty(global_names{index_in_p});
                        if already
                            already = (isempty(strfind(global_names{index_in_p},'internal')) | isempty(strfind(global_names{index_in_p},'ans')));
                        end
                    else
                        already = 0;
                    end
                    
                    if ~isempty(index_in_p) & ~already
                        B = reshape(getbasematrix(thevars,vars(ii)),size(thevars,1),size(thevars,2));
                        [ix,jx,kx] = find(B);
                        ix=ix(1);
                        jx=jx(1);
                        global_names{index_in_p}=[W(i).name '(' num2str(ix) ',' num2str(jx) ')'];
                    end
                end
                
            else
                % Case 3
                vars = getvariables(thevars);
                indicies = find(ismember(vars,global_LinearVariables));
                for ii = indicies
                    index_in_p = find(ismember(global_LinearVariables,vars(ii)));
                    
                    if ~isempty(index_in_p)
                        already = ~isempty(global_names{index_in_p});
                        if already
                            already = (isempty(strfind(global_names{index_in_p},'internal')) | isempty(strfind(global_names{index_in_p},'ans')));
                        end
                    else
                        already = 0;
                    end
                    
                    if ~isempty(index_in_p) & ~already
                        global_names{index_in_p}=[W(i).name '(' num2str(ii) ')'];
                    end
                end
            end
            
        elseif is(thevars,'sdpcone')
            % Case 3
            vars = getvariables(thevars);
            indicies = find(ismember(vars,global_LinearVariables));
            for ii = indicies
                index_in_p = find(ismember(global_LinearVariables,vars(ii)));
                if ~isempty(index_in_p)
                    already = ~isempty(global_names{index_in_p});
                    if already
                        already = ~strfind(global_names{index_in_p},'internal');
                    end
                else
                    already = 0;
                end
                
                if ~isempty(index_in_p) & ~already
                    B = reshape(getbasematrix(thevars,vars(ii)),size(thevars,1),size(thevars,2));
                    [ix,jx,kx] = find(B);
                    ix=ix(1);
                    jx=jx(1);
                    global_names{index_in_p}=[W(i).name '(' num2str(ix) ',' num2str(jx) ')'];
                end
            end
            
        else
            % Case 4
            vars = getvariables(thevars);
            indicies = find(ismember(vars,global_LinearVariables));
            
            for i = indicies
                index_in_p = find(ismember(global_LinearVariables,vars(i)));
                if ~isempty(index_in_p) & isempty(global_names{index_in_p})
                    global_names{index_in_p}=['internal(' num2str(vars(i)) ')'];
                end
            end
            
        end
    end
end

for pi = 1:size(pvec,1)
    for pj = 1:size(pvec,2)
        p = pvec(pi,pj);
        if isa(p,'double')
            symb_p = num2str2(p,precision);
        else            
            symb_p = createSymbolicExpression(p,global_LinearVariables,global_names,precision);           
        end        
        symb_pvec{pi,pj} = symb_p;
    end
end

if prod(size(symb_pvec))==1 & nargout==0
    display(symb_pvec{1,1});
    clear symb_pvec
end

function symb_p =  createSymbolicExpression(p,global_LinearVariables,global_names,precision)
LinearVariables = depends(p);
x = recover(LinearVariables);
[exponent_p,ordered_list] = exponents(p,x);
exponent_p = full(exponent_p);
[~,map] = ismember(LinearVariables,global_LinearVariables);
names = {global_names{map}};

% Remove 0 constant
symb_p = '';
if size(ordered_list,1)>0
    nummonoms = size(ordered_list,1);
    if full(getbasematrix(p,0)) ~= 0
        symb_p = num2str2(full(getbasematrix(p,0)),precision);
    end
elseif all(exponent_p(1,:)==0)
    symb_p = num2str2(full(getbasematrix(p,0)),precision);
    exponent_p = exponent_p(2:end,:);
    nummonoms = size(exponent_p,1);
else
    nummonoms = size(exponent_p,1);
end

% Loop through all monomial terms
for i = 1:nummonoms
    coeff = full(getbasematrixwithoutcheck(p,i));
    switch coeff
        case 1
            coeff='+';
        case -1
            coeff = '-';
        otherwise
            if isreal(coeff)
                if coeff >0
                    coeff = ['+' num2str2(coeff,precision)];
                else
                    coeff=[num2str2(coeff,precision)];
                end
            else
                coeff = ['+' '(' num2str2(coeff,precision) ')' ];
            end
    end
    if isempty(ordered_list)
        symb_p = [symb_p coeff symbmonom(names,exponent_p(i,:))];
    else
        symb_p = [symb_p coeff symbmonom_noncommuting(names,ordered_list(i,:))];
    end
end
% Clean up some left overs, lazy coding...
symb_p = strrep(symb_p,'+*','+');
symb_p = strrep(symb_p,'-*','-');
if symb_p(1)=='+'
    symb_p = symb_p(2:end);
end
if symb_p(1)=='*'
    symb_p = symb_p(2:end);
end


function s = symbmonom(names,monom)
s = '';
for j = 1:length(monom)
    if abs( monom(j))>0
        if isempty(names{j})
            names{j} = ['internal(' num2str(j) ')'];
        end
        s = [s '*' names{j}];
        if monom(j)~=1
            s = [s '^' num2str(monom(j))];
        end
    end
end

function s = symbmonom_noncommuting(names,monom)
s = '';
j = 1;
while j <= length(monom)
    if abs( monom(j))>0
        if isempty(names{monom(j)})
            names{monom(j)} = ['internal(' num2str(j) ')'];
        end
        s = [s '*' names{monom(j)}];
        power = 1; 
        k = j;
        while j<length(monom) & monom(j) == monom(j+1)
            power = power + 1;
            j = j + 1;
            %if j == (length(monom)-1)
            %    j = 5;
            %end
        end
        if power~=1
            s = [s '^' num2str(power)];
        end
    end
    j = j + 1;
end

function s = num2str2(x,precision)
if isinf(precision) 
    s = sprintf('%.12g',x);
else  
    s = sprintf('%.12g',round(x*10^precision)/10^precision);
end
s(s==10)=[];
s(s==32)=[];


