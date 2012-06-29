function Freal = imag2reallmi(F)

% Author Johan Löfberg
% $Id: imag2reallmi.m,v 1.6 2006-12-18 14:42:28 joloef Exp $

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
                               
                if hasfactors(Freal.clauses{j}.data)                    
                    % We use factor-tracking lifting
                   Freal.clauses{j}.data = lift2real(F.clauses{j}.data);                    
                else
                    % Fast without tracking
                    reF=real(F.clauses{j}.data);
                    imF=imag(F.clauses{j}.data);                                            
                    Freal.clauses{j}.data = kron(eye(2),reF) + kron([0 1;-1 0],imF);
                end
                
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