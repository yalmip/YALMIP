function Z = uminus(P)
%display           Overloaded

Z = P;
Z.cx =  - Z.cx;
Z.gain = -Z.gain;
