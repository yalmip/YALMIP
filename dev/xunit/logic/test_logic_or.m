function test_logic_or
binvar y1 y2 y3
solvesdp([true(or(y1,y2,y3))],y1+y2+y3)
mbg_asserttrue(sum(double([y1 y2 y3]))==1);
