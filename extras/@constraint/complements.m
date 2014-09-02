function F = complements(C1,C2)
%COMPLEMENTS Defines complementary constraints
%   
%   F = COMPLEMENTS(C1,C2)   

F = complements(lmi(C1),lmi(C2));
	