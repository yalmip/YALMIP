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

allnames = {};
for pi = 1:size(pvec,1)
    for pj = 1:size(pvec,2)
        Y.type = '()';
        Y.subs = [{pi} {pj}];
        p = subsref(pvec,Y);
      %  p = pvec(pi,pj);
        
        if isa(p,'double')
            symb_p = num2str(p);
        else
            LinearVariables = depends(p);
            x = recover(LinearVariables);
            exponent_p = full(exponents(p,x));
            names = cell(length(x),1);
            for i = 1:length(names)
                names{i} = ['x' num2str(LinearVariables(i))];    
                allnames{end+1} = names{i};
            end
            
            symb_p = '';
            if all(exponent_p(1,:)==0)
                symb_p = num2str(full(getbasematrix(p,0)));
                exponent_p = exponent_p(2:end,:);
            end
           
            for i = 1:size(exponent_p,1)
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
                symb_p = [symb_p coeff symbmonom(names,exponent_p(i,:))];                                
            end
            if symb_p(1)=='+'
                symb_p = symb_p(2:end);
            end
        end
        symb_pvec{pi,pj} = symb_p;
    end
end
allnames = unique(allnames);
for i = 1:length(allnames)
    evalin('caller',['syms ' allnames{i}]);
end


S = '';
for pi = 1:size(pvec,1)
    ss = '';
    for pj = 1:size(pvec,2)
        ss = [ss ' ' symb_pvec{pi,pj} ','];
    end
    S = [S ss ';'];
end
S = ['[' S ']']   ;            
symb_pvec = evalin('caller',S);


function s = symbmonom(names,monom)
s = '';
for j = 1:length(monom)
    if abs( monom(j))>0
        s = [s names{j}];
        if monom(j)~=1
            s = [s '^' num2str(monom(j))];
        end
        s =[s '*'];
    end
    
end
if isequal(s(end),'*')
    s = s(1:end-1);
end

function s = num2str2(x)
s = num2str(full(x));
if isequal(s,'1')
    s = '';
end
if isequal(s,'-1')
    s = '-';
end

        