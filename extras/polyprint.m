function symb_pvec = polyprint(pvec)
%POLYPRINT Pretty print polynomial expression
%
% POLYPRINT is obsolete. Use SDISPLAY instead.

for pi = 1:size(pvec,1)
    for pj = 1:size(pvec,2)
        p = pvec(pi,pj);
        
        if isa(p,'double')
            symb_p = num2str(p);
        else
            LinearVariables = depends(p);
            x = recover(LinearVariables);
            exponent_p = full(exponents(p,x));
            names = cell(length(x),1);
            W = evalin('caller','whos');
            for i = 1:size(W,1)
                if strcmp(W(i).class,'sdpvar')% | strcmp(W(i).class,'lmi')
                    thevars = evalin('caller',W(i).name)    ;
                    if is(thevars,'scalar') & is(thevars,'linear') & length(getvariables(thevars))==1
                        index_in_p = find(ismember(LinearVariables,getvariables(thevars)));
                        if ~isempty(index_in_p)
                            names{index_in_p}=W(i).name;
                        end
                    end
                end
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

function s = symbmonom(names,monom)
s = '';
for j = 1:length(monom)
    if abs( monom(j))>0
        s = [s names{j}];
        if monom(j)~=1
            s = [s '^' num2str(monom(j))];
        end
    end
end

function s = num2str2(x)
        s = num2str(full(x));
        if isequal(s,'1')
            s = '';
        end
        if isequal(s,'-1')
            s = '-';
        end
  
        