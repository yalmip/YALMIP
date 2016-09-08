function [Z,coeffs] = createDualParameterization(UncertaintySet,v,degree);

if nargin == 1
    degree = repmat(2,1,length(UncertaintySet));
elseif length(degree)==1
    degree = repmat(degree,1,length(UncertaintySet));
end
coeffs = [];

for i = 1:length(UncertaintySet)
     n = length(sdpvar(UncertaintySet(i)));
     if is(UncertaintySet(i),'elementwise')    
        [Z{i},coeffsi] = matrixpolynomial(v,[n 1],degree(i));
     elseif is(UncertaintySet(i),'equality')  
        [Z{i},coeffsi] = matrixpolynomial(v,[n 1],degree(i));
     elseif is(UncertaintySet(i),'lmi')      
        [Z{i},coeffsi] = matrixpolynomial(v,[n n],degree(i));
    elseif is(UncertaintySet(i),'socp')
        [Z{i},coeffsi] = matrixpolynomial(v,[n 1],degree(i));   
    end   
    coeffs = [coeffs;coeffsi(:)];
end
coeffs = recover(getvariables(coeffs));