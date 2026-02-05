function S = cov(w)    
% FIXME: Lazy for now, do this explicitly
Ew = expected(w);
Eww = expected(w*w');
S = Eww-Ew*Ew';
