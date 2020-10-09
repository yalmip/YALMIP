function p = propagate_bounds_from_equalities(p)

if p.K.f == 0
    return
end

try
    if ~any(abs(p.ub(p.branch_variables)-p.lb(p.branch_variables))>p.options.bmibnb.vartol)
        return
    end
catch
end

LU = [p.lb p.ub];

p_F_struc = p.F_struc;
n_p_F_struc_cols = size(p_F_struc,2);

fixedVars = find(p.lb == p.ub & p.variabletype(:) == 0);
if ~isempty(fixedVars)
    p_F_struc_forbilin = p_F_struc;
    p_F_struc_forbilin(:,1) = p_F_struc(:,1) + p_F_struc(:,1+fixedVars)*p.lb(fixedVars);
    p_F_struc_forbilin(:,1+fixedVars) = 0;
else
    p_F_struc_forbilin=p_F_struc;
end

usedVariables = find(any(p.F_struc(:,2:end)));
if all(isinf(p.lb(usedVariables))) & all(isinf(p.ub(usedVariables)))
    return
end

if p.K.f >0
    interestingRows = find(p_F_struc(1:p.K.f,1));    
    if ~isempty(interestingRows)
        S = p_F_struc(interestingRows,:);
        S(:,1)=0;
        S = sum(S|S,2) - abs(sum(S,2));
        interestingRows = interestingRows(find(S==0));
        for j = interestingRows(:)'
            thisrow = p_F_struc(j,:);
            if thisrow(1)<0
                thisrow = -thisrow;
            end
            [row,col,val] = find(thisrow);
            % Find bounds from sum(xi) = 1, xi>0
            if all(val(2:end) < 0)
                usedVars = col(2:end)-1;
                if all(p.lb(usedVars)>=0)
                    p.ub(usedVars) = min( p.ub(usedVars) , val(1)./abs(val(2:end)'));
                end
            end
        end
    end
    
    % Presolve from bilinear x*y == k
    if any(p.variabletype == 1)
        for j = 1:p.K.f
            if p_F_struc_forbilin(j,1)~=0
                [row,col,val] = find(p_F_struc_forbilin(j,:));
                
                % Find bounds from sum(xi) = 1, xi>0
                if length(col)==2
                    val = val/val(2); % val(1) + x==0
                    var = col(2)-1;
                    if p.variabletype(var)==1
                        [ij] = find(p.monomtable(var,:));
                        if p.lb(ij(1))>=0 & p.lb(ij(2))>=0
                            % xi*xj == val(1)
                            if -val(1)<0
                                p.feasible = 0;
                                return
                            else
                                p.ub(ij(2)) = min( p.ub(ij(2)),-val(1)/p.lb(ij(1)));
                                p.ub(ij(1)) = min( p.ub(ij(1)),-val(1)/p.lb(ij(2)));
                                p.lb(ij(2)) = max( p.lb(ij(2)),-val(1)/p.ub(ij(1)));
                                p.lb(ij(1)) = max( p.lb(ij(1)),-val(1)/p.ub(ij(2)));
                            end
                        elseif p.ub(ij(1))<=0 & p.ub(ij(2))<=0
                            if -val(1)<0
                                p.feasible = 0;
                                return
                            else
                                x = ij(1);y = ij(2);
                                if p.ub(y)
                                    p.lb(x) = max(p.lb(x),-val(1)/p.ub(y));
                                end
                                if p.ub(x)
                                    p.lb(y) = max(p.lb(y),-val(1)/p.ub(x));
                                end
                                if p.lb(y)
                                    p.ub(x) = min( p.ub(x),-val(1)/p.lb(y));
                                end
                                if p.ub(x)
                                    p.ub(y) = min( p.ub(y),-val(1)/p.lb(x));                                
                                end
                            end
                            
                        elseif -val(1)>0 &  p.lb(ij(1))>=0
                            p.lb(ij(2)) = max(0,p.lb(ij(2)));
                        elseif -val(1)>0 &  p.lb(ij(2))>=0
                            p.lb(ij(1)) = max(0,p.lb(ij(1)));
                        elseif val(1) < 0 &&  val(2) > 0 % x*y == pos so sign can be derived at least
                            if p.lb(ij(1))>=0
                                p.lb(ij(2)) = max(0,p.lb(ij(2)));
                            elseif p.lb(ij(2))>=0
                                p.lb(ij(1)) = max(0,p.lb(ij(1)));
                            elseif p.ub(ij(1))<=0
                                p.ub(ij(2)) = min(0,p.ub(ij(2)));
                            elseif p.ub(ij(2))<=0
                                 p.ub(ij(1)) = min(0,p.ub(ij(1)));
                            end
                        end
                    end
                end
            end
        end
    end
    
   
    A = p.F_struc(1:p.K.f,2:end);
    AT = A';
    Ap = max(0,A);ApT = Ap';
    Am = min(0,A);AmT = Am';
    
    two_terms = sum(p.F_struc(1:p.K.f,2:end) | p.F_struc(1:p.K.f,2:end),2)==2;
      
    for j = find(sum(p.F_struc(1:p.K.f,2:end) | p.F_struc(1:p.K.f,2:end),2)>1)'
        % Simple x == y
        done = 0;
        b = full(p_F_struc(j,1));
        if b==0 & two_terms(j)
            [row,col,val] = find(p_F_struc(j,:));
            if length(row) == 2
                if val(1) == -val(2)
                    p.lb(col(1)-1) = max(p.lb(col(1)-1),p.lb(col(2)-1));
                    p.lb(col(2)-1) = max(p.lb(col(1)-1),p.lb(col(2)-1));
                    p.ub(col(1)-1) = min(p.ub(col(1)-1),p.ub(col(2)-1));
                    p.ub(col(2)-1) = min(p.ub(col(1)-1),p.ub(col(2)-1));
                    done = 1;
                elseif val(1) == val(2)
                    p.lb(col(1)-1) = max(p.lb(col(1)-1),-p.ub(col(2)-1));
                    p.lb(col(2)-1) = max(-p.ub(col(1)-1),p.lb(col(2)-1));
                    p.ub(col(1)-1) = min(p.ub(col(1)-1),-p.lb(col(2)-1));
                    p.ub(col(2)-1) = min(-p.lb(col(1)-1),p.ub(col(2)-1));
                    done = 1;
                end
            end
        end
        if ~done
            a = AT(:,j)';
            ap = (ApT(:,j)');
            am = (AmT(:,j)');
            find_a = find(a);
            
            p_ub = p.ub(find_a);
            p_lb = p.lb(find_a);
            
            inflb = isinf(p_lb);
            infub = isinf(p_ub);
            if ~all(inflb & infub)
                if  any(inflb) | any(infub)
                    [p_lb,p_ub] = propagatewINFreduced(full(a(find_a)),full(ap(find_a)),full(am(find_a)),p_lb,p_ub,b);
                    p.lb(find_a) = p_lb;
                    p.ub(find_a) = p_ub;
                else
                    [p_lb,p_ub] = propagatewoINFreduced(full(a(find_a)),full(ap(find_a)),full(am(find_a)),p_lb,p_ub,b);
                    p.lb(find_a) = p_lb;
                    p.ub(find_a) = p_ub;
                end
            end
        end
    end
end
close = find(abs(p.lb - p.ub) < 1e-12);
p.lb(close) = (p.lb(close)+p.ub(close))/2;
p.ub(close) = p.lb(close);
% Numerical issues easily propagates, so widen weird close to feasible box 
p = widenSmashedBox(p);
p = update_integer_bounds(p);

if ~isequal(LU,[p.lb p.ub])
    p.changedbounds = 1;
end

if any(p.lb > p.ub)
    p.feasible = 0;
end



function [p_lb,p_ub] = propagatewINFreduced(a,ap,am,p_lb,p_ub,b);
%a = AT(:,j)';
%ap = (ApT(:,j)');
%am = (AmT(:,j)');

%p_ub = p.ub;
%p_lb = p.lb;

%find_a = find(a);
%  find_a = find_a(min(find(isinf(p.lb(find_a)) | isinf(p.ub(find_a)))):end);
for k = 1:length(a)%find_a
    
    p_ub_k = p_ub(k);
    p_lb_k = p_lb(k);
    
    if (p_ub_k-p_lb_k) > 1e-8
        L = p_lb;
        U = p_ub;
        L(k) = 0;
        U(k) = 0;
        ak = a(k);
        
        if ak < 0
            ak = -ak;
            aa = am;
            am = -ap;
            ap = -aa;
            b = -b;
            a = -a;
        end
        
        if ak > 0       
            use1 = find(ap'~=0);
            use2 = find(am'~=0);            
            newlower = (-b - ap(use1)*U(use1) - am(use2)*L(use2))/ak;
            newupper = (-b - am(use2)*U(use2) - ap(use1)*L(use1))/ak;
            %newlower = (-b - ap*U - am*L)/ak;
            %newupper = (-b - am*U - ap*L)/ak;
        else
            newlower = (-b - am*U - ap*L)/ak;
            newupper = (-b - ap*U - am*L)/ak;
        end       
        if p_ub_k>newupper
            p_ub(k) = newupper;           
        end
        if p_lb_k<newlower
            p_lb(k) = newlower;
        end
    end
end
p.ub = p_ub;
p.lb = p_lb;

function [p_lb,p_ub] = propagatewoINFreduced(a,ap,am,p_lb,p_ub,b);

L = p_lb;
U = p_ub;

apU = ap*U;
amU = am*U;
apL = ap*L;
amL = am*L;

papU = ap.*U';
pamU = am.*U';
papL = ap.*L';
pamL = am.*L';

minusbminusapUminusamL = -b-apU-amL;
minusbminusamUminusapL = -b-amU-apL;
for k = 1:length(a)%find_a
    
    p_ub_k = p_ub(k);
    p_lb_k = p_lb(k);
    
    if (p_ub_k-p_lb_k) > 1e-8
        ak = a(k);
        if ak > 0
            %newlower = (-b-apU+papU(k)-amL+pamL(k) )/ak;
            %newupper = (-b-amU+pamU(k)-apL+papL(k) )/ak;
            newlower = -1e-15 + (minusbminusapUminusamL+papU(k)+pamL(k) )/ak;
            newupper = 1e-15 + (minusbminusamUminusapL+pamU(k)+papL(k) )/ak;
        else
            newlower = -1e-15 + (minusbminusamUminusapL+pamU(k)+papL(k) )/ak;
            newupper = 1e-15 + (minusbminusapUminusamL+papU(k)+pamL(k) )/ak;
        end
        if p_ub_k>newupper
            p_ub(k) = newupper;
            U(k) = newupper;
            apU = ap*U;
            amU = am*U;
            papU = ap.*U';
            pamU = am.*U';
            minusbminusapUminusamL = -b-apU-amL;
            minusbminusamUminusapL = -b-amU-apL;
        end
        if p_lb_k<newlower
            p_lb(k) = newlower;
            L(k) = newlower;
            apL = ap*L;
            amL = am*L;
            papL = ap.*L';
            pamL = am.*L';
            minusbminusapUminusamL = -b-apU-amL;
            minusbminusamUminusapL = -b-amU-apL;
        end
    end
end
%p.ub = p_ub;
%p.lb = p_lb;

function p = propagatewINF(p,AT,ApT,AmT,j,b);
a = AT(:,j)';
ap = (ApT(:,j)');
am = (AmT(:,j)');


p_ub = p.ub;
p_lb = p.lb;

find_a = find(a);
%  find_a = find_a(min(find(isinf(p.lb(find_a)) | isinf(p.ub(find_a)))):end);
for k = find_a
    
    p_ub_k = p_ub(k);
    p_lb_k = p_lb(k);
    
    if (p_ub_k-p_lb_k) > 1e-8
        L = p_lb;
        U = p_ub;
        L(k) = 0;
        U(k) = 0;
        ak = a(k);
        if ak > 0
            newlower = (-b - ap*U - am*L)/ak;
            newupper = (-b - am*U - ap*L)/ak;
        else
            newlower = (-b - am*U - ap*L)/ak;
            newupper = (-b - ap*U - am*L)/ak;
        end
        %                     if isinf(newlower) | isinf(newupper)
        %                         z = newlower;
        %                     end
        if p_ub_k>newupper
            p_ub(k) = newupper;           
        end
        if p_lb_k<newlower
            p_lb(k) = newlower;
        end
    end
end
p.ub = p_ub;
p.lb = p_lb;

function p = propagatewoINF(p,AT,ApT,AmT,j,b);
a = full(AT(:,j)');
ap = full((ApT(:,j)'));
am = full((AmT(:,j)'));


p_ub = p.ub;
p_lb = p.lb;

find_a = find(a);

    L = p_lb;
    U = p_ub;
    
    apU = ap*U;
    amU = am*U;
    apL = ap*L;
    amL = am*L;
    
    papU = ap.*U';
    pamU = am.*U';
    papL = ap.*L';
    pamL = am.*L';
    
    
    minusbminusapUminusamL = -b-apU-amL;    
    minusbminusamUminusapL = -b-amU-apL;    
for k = find_a
    
    p_ub_k = p_ub(k);
    p_lb_k = p_lb(k);
    
    if (p_ub_k-p_lb_k) > 1e-8
        ak = a(k);
        if ak > 0
            %newlower = (-b-apU+papU(k)-amL+pamL(k) )/ak;
            %newupper = (-b-amU+pamU(k)-apL+papL(k) )/ak;
            newlower = (minusbminusapUminusamL+papU(k)+pamL(k) )/ak;
            newupper = (minusbminusamUminusapL+pamU(k)+papL(k) )/ak;
        else
            newlower = (minusbminusamUminusapL+pamU(k)+papL(k) )/ak;
            newupper = (minusbminusapUminusamL+papU(k)+pamL(k) )/ak;
        end
        if p_ub_k>newupper
            p_ub(k) = newupper;
            U(k) = newupper;
            apU = ap*U;
            amU = am*U;
            papU = ap.*U';
            pamU = am.*U';
            minusbminusapUminusamL = -b-apU-amL;    
            minusbminusamUminusapL = -b-amU-apL;    
        end
        if p_lb_k<newlower
            p_lb(k) = newlower;
            L(k) = newlower;
            apL = ap*L;
            amL = am*L;
            papL = ap.*L';
            pamL = am.*L';
            minusbminusapUminusamL = -b-apU-amL;    
            minusbminusamUminusapL = -b-amU-apL;    
        end
    end
end
p.ub = p_ub;
p.lb = p_lb;






























