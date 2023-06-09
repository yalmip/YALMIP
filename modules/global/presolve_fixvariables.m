function [fix_up,fix_down] = presolve_fixvariables(A,b,c,lb,ub,monotinicity)

% These are optimally (or w.l.o.g) set to upper bound
not_in_obj = find(c<=0);
% Setting to 1 makes Ax<b increase feasible set for these variables
constrained_blow = all(-A(:,not_in_obj)>=0,1);
% and they enter via a psd matrix in all sdp constraints
sdp_positive = monotinicity(not_in_obj) == -1;
% these variables satisffy all constraints
can_fix = not_in_obj(find(constrained_blow & sdp_positive));
% these variables are still not fixed
still_on = find(lb==0 & ub==1);
% so we can fix these
fix_up = intersect(can_fix,still_on);

% These are optimally (or w.l.o.g) set to lower bound
not_in_obj = find(c>=0);
% Setting to 1 makes Ax<b increase feasible set for these variables
constrained_blow = all(A(:,not_in_obj)>=0,1);
% and they enter via a psd matrix in all sdp constraints
sdp_positive = monotinicity(not_in_obj) == 1;
% these variables satisffy all constraints
can_fix = not_in_obj(find(constrained_blow & sdp_positive));
% these variables are still not fixed
still_on = find(lb==0 & ub==1);
% so we can fix these
fix_down = intersect(can_fix,still_on);

