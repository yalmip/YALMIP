function y = sqrtm_internal(x)

if x>=0
    y = sqrt(x);
else
    y = -x.^2;
end
