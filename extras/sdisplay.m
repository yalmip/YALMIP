function symb_pvec = sdisplay(pvec,symbolicname)
%SDISPLAY Symbolic display of SDPVAR expression
%
% Note that the symbolic display only work if all
% involved variables are explicitely defined as
% scalar variables.
%
% Variables that not are defined as scalars
% will be given the name ryv(i). ryv means
% recovered YALMIP variables, i indicates the
% index in YALMIP (i.e. the result from getvariables)
%
% If you want to change the generic name ryv, just
% pass a second string argument
%
% EXAMPLES
%  sdpvar x y
%  sdisplay(x^2+y^2)
%    ans =
%       'x^2+y^2'
%
%  t = sdpvar(2,1);
%  sdisplay(x^2+y^2+t'*t)
%    ans =
%      'x^2+y^2+ryv(5)^2+ryv(6)^2'

r1=1:size(pvec,1);
r2=1:size(pvec,2);

for pi = 1:size(pvec,1)
    for pj = 1:size(pvec,2)

        p = pvec(pi,pj);

        if isa(p,'double')
            symb_p = num2str(p,10);
        else
            LinearVariables = depends(p);
            x = recover(LinearVariables);
            [exponent_p,ordered_list] = exponents(p,x);
            exponent_p = full(exponent_p);
            names = cell(length(x),1);

            % First, some boooring stuff. we need to
            % figure out the symbolic names and connect
            % these names to YALMIPs variable indicies
            W = evalin('caller','whos');
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
                        index_in_p = find(ismember(LinearVariables,getvariables(thevars)));

                        if ~isempty(index_in_p)
                            already = ~isempty(names{index_in_p});
                            if already
                                already = ~strfind(names{index_in_p},'internal');
                                if isempty(already)
                                    already = 0;
                                end
                            end
                        else
                            already = 0;
                        end

                        if ~isempty(index_in_p) & ~already
                            % Case 1
                            names{index_in_p}=W(i).name;
                        end
                    elseif is(thevars,'lpcone')

                        if size(thevars,1)==size(thevars,2)
                            % Case 2
                            vars = getvariables(thevars);
                            indicies = find(ismember(vars,LinearVariables));
                            for ii = indicies
                                index_in_p = find(ismember(LinearVariables,vars(ii)));

                                if ~isempty(index_in_p)
                                    already = ~isempty(names{index_in_p});
                                    if already
                                        already = ~strfind(names{index_in_p},'internal');
                                        if isempty(already)
                                    already = 0;
                                end
                                    end
                                else
                                    already = 0;
                                end

                                if ~isempty(index_in_p) & ~already
                                    B = reshape(getbasematrix(thevars,vars(ii)),size(thevars,1),size(thevars,2));
                                    [ix,jx,kx] = find(B);
                                    ix=ix(1);
                                    jx=jx(1);
                                    names{index_in_p}=[W(i).name '(' num2str(ix) ',' num2str(jx) ')'];
                                end
                            end

                        else
                            % Case 3
                            vars = getvariables(thevars);
                            indicies = find(ismember(vars,LinearVariables));
                            for ii = indicies
                                index_in_p = find(ismember(LinearVariables,vars(ii)));

                                if ~isempty(index_in_p)
                                    already = ~isempty(names{index_in_p});
                                    if already
                                        already = ~strfind(names{index_in_p},'internal');
                                        if isempty(already)
                                    already = 0;
                                end
                                    end
                                else
                                    already = 0;
                                end

                                if ~isempty(index_in_p) & ~already
                                    names{index_in_p}=[W(i).name '(' num2str(ii) ')'];
                                end
                            end
                        end

                    elseif is(thevars,'sdpcone')
                        % Case 3
                        vars = getvariables(thevars);
                        indicies = find(ismember(vars,LinearVariables));
                        for ii = indicies
                            index_in_p = find(ismember(LinearVariables,vars(ii)));
                            if ~isempty(index_in_p)
                                already = ~isempty(names{index_in_p});
                                if already
                                    already = ~strfind(names{index_in_p},'internal');
                                end
                            else
                                already = 0;
                            end

                            if ~isempty(index_in_p) & ~already
                                B = reshape(getbasematrix(thevars,vars(ii)),size(thevars,1),size(thevars,2));
                                [ix,jx,kx] = find(B);
                                ix=ix(1);
                                jx=jx(1);
                                names{index_in_p}=[W(i).name '(' num2str(ix) ',' num2str(jx) ')'];
                            end
                        end

                    else
                        % Case 4
                        vars = getvariables(thevars);
                        indicies = find(ismember(vars,LinearVariables));

                        for i = indicies
                            index_in_p = find(ismember(LinearVariables,vars(i)));
                            if ~isempty(index_in_p) & isempty(names{index_in_p})
                                names{index_in_p}=['internal(' num2str(vars(i)) ')'];
                            end
                        end

                    end
                end
            end

            % Okay, now got all the symbolic names compiled.
            % Time to construct the expression

            % The code below is also a bit fucked up at the moment, due to
            % the experimental code with noncommuting stuff 
            
            % Remove 0 constant
            symb_p = '';
            if size(ordered_list,1)>0
                nummonoms = size(ordered_list,1);
                if full(getbasematrix(p,0)) ~= 0
                     symb_p = num2str(full(getbasematrix(p,0)));
                end
            elseif all(exponent_p(1,:)==0)
                symb_p = num2str(full(getbasematrix(p,0)),10);
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
                                coeff = ['+' num2str2(coeff)];
                            else
                                coeff=[num2str2(coeff)];
                            end
                        else
                            coeff = ['+' '(' num2str2(coeff) ')' ];
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
        end
        symb_pvec{pi,pj} = symb_p;
    end
end

if prod(size(symb_pvec))==1 & nargout==0
    display(symb_pvec{1,1});
    clear symb_pvec
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

function s = num2str2(x)
s = evalc('disp(x)');
s(s==10)=[];
s(s==32)=[];
%disp(full(x))
%s = num2str(full(x),10);
if isequal(s,'1')
    s = '';
end
if isequal(s,'-1')
    s = '-';
end

