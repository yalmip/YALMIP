function mono = decreasing_except_at(xL,xU,p)

if xU <= p || xL >= p
    mono = 'decreasing';
else
    mono = 'none';
end