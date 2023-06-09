function mono = increasing_except_at(xL,xU,p)

if xU <= p || xL >= p
    mono = 'increasing';
else
    mono = 'none';
end