function reDeclareForBinaryMax(y,X);
vX = getvariables(X);
B = getbase(X);
vB = yalmip('binvariables');
if all(ismember(vX,vB))
    if all(B(:,1)==0)
        if is(X,'lpcone') | is(X,'sdpcone')
            yalmip('setbinvariables',[vB getvariables(y)]);
        end
    end
end