function symb_pvec = sdisplay2(p,names)

% TODO: sdpvar a b
%       x = (a+b)^0.3        -- writes "mpower_internal"
%
% TODO: sdpvar h x k
%       sdisplay2(h*x - k)   -- misses the minus in front of "k"

if nargin < 2
    LinearVariables = 1:yalmip('nvars');
    x = recover(LinearVariables);
    names = cell(length(x),1);

    W = evalin('caller','whos');
    for i = 1:size(W,1)
        if strcmp(W(i).class,'sdpvar') | strcmp(W(i).class,'ncvar')
            % Get the SDPVAR variable
            thevars = evalin('caller',W(i).name);

            % Distinguish 4 cases
            % 1: Sclalar varible x
            % 2: Vector variable x(i)
            % 3: Matrix variable x(i,j)
            % 4: Variable not really defined
            if is(thevars,'scalar') & is(thevars,'linear') & length(getvariables(thevars))==1 & isequal(getbase(thevars),[0 1])
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
end

[mt,vt] = yalmip('monomtable');
ev = yalmip('extvariables');

for i = 1:size(p, 1)
    for j = 1:size(p, 2)
        symb_pvec{i, j} = symbolicdisplay(p(i, j), names, vt, ev, mt);
    end
end


%------------------------------------------------------------------------
function expression = symbolicdisplay(p,names,vt,ev,mt)

sp = size(p);
if any(sp > 1)
    out = '[';
else
    out = '';
end
p_orig = p;
for i1 = 1:sp(1)
    for i2 = 1:sp(2)
        p = p_orig(i1, i2);
        basis = getbase(p);
        if basis(1)~=0
            expression = [num2str(basis(1)) '+'];
        else
            expression = [''];
        end

        [dummy, variables, coeffs] = find(basis(2:end));
        variables = getvariables(p);
        for i = 1:length(coeffs)
            if coeffs(i)==1
                expression = [expression symbolicmonomial(variables(i), ...
                    names,vt,ev,mt) '+'];
            else
                expression = [expression num2str(coeffs(i)) '*' ...
                    symbolicmonomial(variables(i),names,vt,ev,mt) '+'];
            end
        end
        expression(end) = [];
        out = [out expression ','];
    end
    out(end) = ';';
end
out(end) = [];
if any(sp > 1)
    out = [out ']'];
end
expression = out;

%------------------------------------------------------------------------
function s = symbolicmonomial(variable,names,vt,ev,mt)

terms = find(mt(variable,:));
if ismember(variable,ev)
    q = yalmip('extstruct',variable);
    s = [q.fcn '(' symbolicdisplay(q.arg{1},names,vt,ev,mt)];
    for i = 2:length(q.arg)-1
        s = [s ',' symbolicdisplay(q.arg{i}, names, vt, ev, mt)];
    end
    s = [s ')'];

elseif ~vt(variable)
    % Linear expression
    s = names{variable};
else
    % Fancy display of a monomial
    s = [''];
    for i = 1:length(terms)
        if mt(variable,terms(i)) == 1
            exponent = '';
        else
            exponent = ['^' num2str(mt(variable,terms(i)))];
        end
        s = [s symbolicmonomial(terms(i),names,vt,ev,mt) exponent '*'];
    end
    s(end)=[];
end
% s = strrep(s,'^1+','+');
% s = strrep(s,'^1*','*');
