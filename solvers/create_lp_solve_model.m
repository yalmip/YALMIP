function lp = create_lp_solve_model(A,b,f,xint,LB,UB,e,options);

[m,n] = size(A);
lp = lp_solve('make_lp', m, n);

lp_solve('set_mat', lp, A);
lp_solve('set_rh_vec', lp, b);
lp_solve('set_obj_fn', lp, f);
lp_solve('set_maxim', lp); % default is solving minimum lp.
for i = 1:length(e)
    if e(i) < 0
        con_type = 1;
    elseif e(i) == 0
        con_type = 3;
    else
        con_type = 2;
    end
    lp_solve('set_constr_type', lp, i, con_type);
end
for i = 1:length(LB)
    %if ~isinf(LB(i))
        lp_solve('set_lowbo', lp, i, LB(i));    
    %end
end
for i = 1:length(UB)    
   if ~isinf(UB(i))
        lp_solve('set_upbo', lp, i, UB(i));    
   end
end
for i = 1:length(xint)
    lp_solve('set_int', lp, xint(i), 1);
end

if options.lpsolve.scalemode~=0
    lp_solve('set_scaling', lp, scalemode);
end

% for i = 1:length(sos)
%     lp_solve('add_SOS', lp, ['dummy' num2str(i)], 1, i, sos{i}, 1:length(sos{i}));
% end


switch options.verbose
    case 0
        lp_solve('set_verbose', lp, 0);%options.verbose)
    case 1
        lp_solve('set_verbose', lp, 4);%options.verbose)
    case 2
        lp_solve('set_verbose', lp, 5);%options.verbose)
    otherwise
        lp_solve('set_verbose', lp, 6);%options.verbose)
end
