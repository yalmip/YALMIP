function F = subsref(F,Y)
%subsref           Overloaded indexing

F = subsref(lmi(F),Y);