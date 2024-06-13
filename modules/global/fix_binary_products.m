function x = set_binary_products(p,x)
if ~isempty(p.binaryProduct)
    x(p.binaryProduct(:,1)) = prod(x(p.binaryProduct(:,2:3)),2);
end