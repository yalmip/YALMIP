function monom = recovermonoms(newton_m,x);

if isempty(newton_m)
    monom = 1;
else
    error('Report this bug (call to recovermonoms with double)');
end