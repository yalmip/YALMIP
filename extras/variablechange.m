function [exponent_p_monoms_changed,varchange] = variablechange(exponent_p_monoms)
%VARIABLECHANGE Internal function to reduce monomials in SOS programs

n = size(exponent_p_monoms,2);
exponent_p_monoms_changed = exponent_p_monoms;
varchange = ones(1,n);

% UGLY!!!!!!!!!!!!!!!!!! Loop and clean up...

% Even multiple of 4?
for i = 1:size(exponent_p_monoms,2)
    rems = rem(exponent_p_monoms(:,i),4);
    if all(rems==0) & all(even(exponent_p_monoms(:,i)/4))
        exponent_p_monoms_changed(:,i)=exponent_p_monoms_changed(:,i)/4;
        varchange(i)=4*varchange(i);
    end
end

% Multiple of 3?
for i = 1:size(exponent_p_monoms,2)
    rems = rem(exponent_p_monoms(:,i),3);
    if all(rems==0)
        exponent_p_monoms_changed(:,i)=exponent_p_monoms_changed(:,i)/3;
        varchange(i)=3;
    end       
end

exponent_p_monoms = exponent_p_monoms_changed;

% Even multiple of 2?
for i = 1:size(exponent_p_monoms,2) 
    rems = rem(exponent_p_monoms(:,i),2);
    if all(rems==0) & all(even(exponent_p_monoms(:,i)/2))
        exponent_p_monoms_changed(:,i)=exponent_p_monoms_changed(:,i)/2;
        varchange(i)=2*varchange(i);
    end
end

return

n = size(exponent_p,2);
exponent_p_changed = exponent_p;
varchange = ones(1,n);

% UGLY!!!!!!!!!!!!!!!!!! Loop and clean up...

% Even multiple of 4?
for i = 1:length(monom_indicies)    
    rems = rem(exponent_p(:,monom_indicies(i)),4);
    if all(rems==0) & all(even(exponent_p(:,monom_indicies(i))/4))
        exponent_p_changed(:,monom_indicies(i))=exponent_p_changed(:,monom_indicies(i))/4;
        varchange(monom_indicies(i))=4*varchange(monom_indicies(i));
    end
end

% Multiple of 3?
for i = 1:length(monom_indicies)
    rems = rem(exponent_p(:,monom_indicies(i)),3);
    if all(rems==0)
        exponent_p_changed(:,monom_indicies(i))=exponent_p_changed(:,monom_indicies(i))/3;
        varchange(monom_indicies(i))=3;
    end       
end

exponent_p = exponent_p_changed;

% Even multiple of 2?
for i = 1:length(monom_indicies)    
    rems = rem(exponent_p(:,monom_indicies(i)),2);
    if all(rems==0) & all(even(exponent_p(:,monom_indicies(i))/2))
        exponent_p_changed(:,monom_indicies(i))=exponent_p_changed(:,monom_indicies(i))/2;
        varchange(monom_indicies(i))=2*varchange(monom_indicies(i));
    end
end
