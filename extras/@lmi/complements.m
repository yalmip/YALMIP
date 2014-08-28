function sys = complements(C1,C2)
%COMPLEMENTS Defines complementary constraints
%   
%   F = COMPLEMENTS(C1,C2)
% 
% Example: F = complements(x >= 0, A*x <= b)
  
if ~(isa(C1,'lmi') & isa(C2,'lmi'))
    error('both arguments in complements must be linear (in)equalities')
end

C1 = flatten(C1);
C2 = flatten(C2);
sys = C1;
sys.clauses{1}.data = [C1.clauses{1}.data(:) C2.clauses{1}.data(:)];
sys.clauses{1}.extra.indicators = binvar(length(C1.clauses{1}.data(:)),1);
sys.clauses{1}.type = 55;