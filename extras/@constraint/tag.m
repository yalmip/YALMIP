function F = tag(F,t)

if nargin > 1
    F = tag(lmi(F),t);
else
    F = tag(lmi(F));
end

