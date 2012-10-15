function F = lowrank(F,x)
% LOWRANK is used to declare that a semidefinite constraint uses data with low rank.
% Used in combination with the solver SDPLR

if nargin > 1
    F = lowrank(lmi(F),x);
else
    F = lowrank(lmi(F));
end
