function newvariabletype = newvariabletypegen(newmt)
newvariabletype = spalloc(size(newmt,1),1,0)';
nonlinear = ~(sum(newmt,2)==1 & sum(newmt~=0,2)==1);
if ~isempty(nonlinear)   
    newvariabletype(nonlinear) = 3;
    quadratic = sum(newmt,2)==2;
    newvariabletype(quadratic) = 2;
    bilinear = max(newmt,[],2)<=1;
    newvariabletype(bilinear & quadratic) = 1;
    sigmonial = any(0>newmt,2) | any(newmt-fix(newmt),2);
    newvariabletype(sigmonial) = 4;
end