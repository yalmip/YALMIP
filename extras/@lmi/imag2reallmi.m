function Freal = imag2reallmi(F)

F = flatten(F);
Counter = size(F.clauses,2);
Freal = F;
j=1;
for i = 1:Counter
    if isreal(F.clauses{i}.data)
        Freal.clauses{j}=F.clauses{i};
        j = j+1;
    else
        switch F.clauses{i}.type
            case 1                
                reF=real(F.clauses{j}.data);
                imF=imag(F.clauses{j}.data);                                                               
                Freal.clauses{j}.data = [reF imF;-imF reF];                                
                j = j+1;
            case {2,3}
                reF=real(F.clauses{j}.data);
                imF=imag(F.clauses{j}.data);
                Freal.clauses{j}.data=[reF;imF];
                j = j+1;
            case 4
                reF=real(F.clauses{j}.data);
                imF=imag(F.clauses{j}.data);
                Freal.clauses{j}.data=[reF(1);reF(2:end);imF(2:end)];
                j = j+1;
            otherwise
                error('Internal bug. Please report.');
        end
    end
end