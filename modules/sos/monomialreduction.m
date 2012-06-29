function exponent_m = monomialreduction(exponent_m,exponent_p,options,csclasses,LPmodel)
%MONOMIALREDUCTION  Internal function for monomial reduction in SOS programs

% Author Johan Löfberg
% $Id: monomialreduction.m,v 1.2 2006-09-26 14:28:43 joloef Exp $


% **********************************************
% TRIVIAL REDUCTIONS (stupid initial generation)
% **********************************************
mindegrees = min(exponent_p,[],1);
maxdegrees = max(exponent_p,[],1);
mindeg = min(sum(exponent_p,2));
maxdeg = max(sum(exponent_p,2));
if size(exponent_m{1},2)==0 % Stupid case : set(sos(parametric))
   if options.verbose>0;disp('Initially 1 monomials in R^0');end
else
    if options.verbose>0;disp(['Initially ' num2str(sum(cellfun('prodofsize',exponent_m)/size(exponent_m{1},2))) ' monomials in R^' num2str(size(exponent_p,2))]);end
end
    
for i = 1:length(csclasses)  
    t = cputime;
    % THE CODE BELOW IS MESSY TO HANDLE SEVERAL BUGS IN MATLAB
    %too_large_term = any(exponent_m-repmat(maxdegrees/2,size(exponent_m,1),1)>0,2);% DOES NOT HANDLE ODD
    % POLYNIMIALS CORRECTLY
    a1 = full(ceil((1+maxdegrees)/2)); % 6.5.1 in linux freaks on sparse stuff...
    if isempty(a1)
        a1 = zeros(size(maxdegrees));
    end
    a2 = full(size(exponent_m{i},1));
    too_large_term = any(exponent_m{i}-repmat(a1,a2,1)>0,2);
    %too_small_term = any(exponent_m-repmat(mindegrees/2,size(exponent_m,1),1)<0,2);
    a1 = full(floor(mindegrees/2));
     if isempty(a1)
        a1 = zeros(size(mindegrees));
    end
    a2 = full(size(exponent_m{i},1));
    too_small_term = any(exponent_m{i}-repmat(a1,a2,1)<0,2);%x^2+xz
    %too_large_sum = any(sum(exponent_m,2)-maxdeg/2>0,2); % DOES NOT HANDLE ODD
    % POLYNIMIALS CORRECTLY
    too_large_sum = any(sum(exponent_m{i},2)-ceil((1+maxdeg)/2)>0,2);
    too_small_sum = any(sum(exponent_m{i},2)-mindeg/2<0,2);
    keep = setdiff1D((1:size(exponent_m{i},1)),find(too_large_term | too_small_term | too_large_sum | too_small_sum));
    exponent_m{i} = exponent_m{i}(keep,:);
    t = cputime-t;    
end
if options.verbose>1;disp(['Removing large/small............Keeping ' num2str(sum(cellfun('prodofsize',exponent_m)/size(exponent_m{1},2))) ' monomials (' num2str(t) 'sec)']);end

% ************************************************
% Homogenuous?
% ************************************************
if all(sum(exponent_p,2)==sum(exponent_p(1,:)))
    for i = 1:length(csclasses)
        j = csclasses{i};
        t = cputime;
        exponent_m{i} = exponent_m{i}(sum(exponent_m{i},2)==sum(exponent_p(1,:))/2,:);
        t = cputime-t;        
    end
    if options.verbose>1;disp(['Homogenuous form!...............Keeping ' num2str(sum(cellfun('prodofsize',exponent_m)/size(exponent_m{1},2))) ' monomials (' num2str(t) 'sec)']);end
end

% ************************************************
% DIAGONAL CONSISTENCY : MONOMIAL ONLY IN
% DIAGONAL, CONSTRAINED TO BE ZER0, CAN BE REMOVED
% ************************************************
if (options.sos.inconsistent==1) &  ~options.sos.csp
    t = cputime;
    keep = consistent(exponent_m{1},exponent_p);
    exponent_m{1} = exponent_m{1}(keep,:);
    t = cputime-t;
    if options.verbose>0;disp(['Diagonal inconsistensies........Keeping ' num2str(size(exponent_m{1},1)) ' monomials (' num2str(t) 'sec)']);end
end

% ***********************************************************
% NEWTON POLYTOPE CHECK
% ***********************************************************
if options.sos.newton
    t = cputime;
    [exponent_m,changed,no_lp_solved] = newtonreduce(exponent_m,exponent_p,options,LPmodel);
    t = cputime-t;
    
    if options.verbose>1 | (options.verbose>0 & 1+changed);
        info = ['Newton polytope (' num2str(no_lp_solved) ' LPs)'];
        info = [info repmat('.',1,32-length(info))];
        disp([info 'Keeping ' num2str(size(exponent_m{end},1)) ' monomials (' num2str(t) 'sec)']);
    end
end

if (options.sos.inconsistent==2) &  ~options.sos.csp
    t = cputime;
    keep = consistent(exponent_m{1},exponent_p);
    exponent_m{1} = exponent_m{1}(keep,:);
    t = cputime-t;
    if options.verbose>1;disp(['Diagonal inconsistensies........Keeping ' num2str(size(exponent_m,1)) ' monomials (' num2str(t) 'sec)']);end
end


