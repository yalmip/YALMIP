function LU = getbounds(F,avoidequalitybounds)

if nargin == 1
    LU = getbounds(lmi(F));
else
    LU = getbounds(lmi(F),avoidequalitybounds);
end
