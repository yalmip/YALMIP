function exponent_m = monomialgeneration(exponent_p,csclasses)
%MONOMIALGENERATION  Internal function to create candidate monomials in SOS programs

exponent_m = [];
for i = 1:length(csclasses)
    if isempty(exponent_p)
        exponent_m{i,1} = zeros(1,0);
    else
        % Create initial set of monomials
        mindegrees = min(exponent_p(:,csclasses{i}));
        maxdegrees = max(exponent_p(:,csclasses{i}));
        if any(2*floor((maxdegrees/2))>maxdegrees)
            error('Highest degree in a variable is odd => not PSD')
        end

        % Make initial generation smarter...        
        exponent_m_temp1 = monolistcoeff(size(csclasses{i},2),ceil(maxdegrees/2),max(ceil(sum(exponent_p,2)/2)));

        [ii,jj] = sort(sum(exponent_m_temp1,2));
        exponent_m_temp1 = exponent_m_temp1(jj,:);
               
        exponent_m_temp2 = zeros(size(exponent_m_temp1,1),size(exponent_p,2));
        exponent_m_temp2(:,csclasses{i}) = exponent_m_temp1;
        exponent_m{i,1} = exponent_m_temp2;
    end
end