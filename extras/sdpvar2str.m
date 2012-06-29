function symb_pvec = sdpvar2str(pvec)
%SDPVAR2STR Converts an SDPVAR object to MATLAB string representation
%
% S = SDPVAR2STR(P)
%
% S : String
% P : SDPVAR object


for pi = 1:size(pvec,1)
    for pj = 1:size(pvec,2)
        p = pvec(pi,pj);
        
        if isa(p,'double')
            symb_p = num2str(p);
        else
            LinearVariables = depends(p);
            x = recover(LinearVariables);
            exponent_p = full(exponents(p,x));
            names = cell(length(LinearVariables),1);
            for i = 1:length(LinearVariables)
                names{i}=['x(' num2str(LinearVariables(i)) ')'];
            end
            
            symb_p = '';
            if all(exponent_p(1,:)==0)
                symb_p = num2str(getbasematrix(p,0));
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
       
        symb_p = strrep(symb_p,'*^0','');
        symb_p = strrep(symb_p,'^0','');
        symb_p = strrep(symb_p,'+*','+');
        symb_p = strrep(symb_p,'-*','-');
        symb_pvec{pi,pj} = symb_p;
    end
end

function s = symbmonom(names,monom)
s = '';
for j = 1:length(monom)
    if monom(j)~=0
        if strcmp(s,'')
            s = [s names{j}];
        else
            s = [s '*' names{j}];
        end
    end
    if monom(j)~=1
        s = [s '^' num2str(monom(j))];
    end
%      if monom(j)>1
%        s = [s '^' num2str(monom(j))];
%    end
    
end

function s = num2str2(x)
        s = num2str(x);
        if isequal(s,'1')
            s = '';
        end
        if isequal(s,'-1')
            s = '-';
        end
  