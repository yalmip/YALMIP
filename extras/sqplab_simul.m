function [outdic,f,ci,ce,cs,g,ai,ae] = sqplab_simul(indic,x,lm,v)

switch indic
    case 0
    case 1
        outdic = 0;
    case 2
        [f,g] = sqplab_fun(x);
        [ci,ce,ai,ae] = sqplab_con(x);
         outdic = 0;
    case 3
        [f,g] = sqplab_fun(x);
        [ci,ce,ai,ae] = sqplab_con(x);
         outdic = 0;
    case 4
        [f,g] = sqplab_fun(x);
        [ci,ce,ai,ae] = sqplab_con(x);
     %   ci = -ci;
        ai = -ai;
        ai = reshape(ai,[],length(g));
        ae = reshape(ae,[],length(g));
        cs = []; 
        outdic = 0;
end
