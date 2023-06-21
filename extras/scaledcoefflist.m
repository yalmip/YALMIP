function c = scaledcoefflist(x,dmax)

    powers = monpowers(length(x),dmax);  
    fact = factorial(powers);
    c = factorial(sum(powers,2)) ./ prod(fact, 2);
    c = arrayfun(@(x) sqrt(x), c);

end