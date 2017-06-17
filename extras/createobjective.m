function [c,Q,f,onlyfeasible] = createobjective(h,G,options,quad_info)
%CREATEOBJECTIVE Internal function to extract data related to objective function

onlyfeasible = 0;
nvars = yalmip('nvars'); 
if isempty(h)
    c=zeros(nvars,1);
    Q = spalloc(nvars,nvars,0);
    f = 0;
    if isempty(G)
        onlyfeasible = 1;
    end
else  
    [n,m]=size(h);
    if n*m>1
        error('Scalar expression to minimize please.');
    else
        % Should check quadratic!!
        if ~(options.relax == 1 | options.relax == 3) & ~isempty(quad_info)
            Qh = quad_info.Q;
            ch = quad_info.c;
            f = quad_info.f;
            xvar = quad_info.x;
            lmi_variables = getvariables(xvar);
            c = zeros(nvars,1);
            Q = spalloc(nvars,nvars,0); 
            c(lmi_variables)=ch;
%             for i=1:length(lmi_variables)
%                 c(lmi_variables(i))=ch(i);
%             end
            if nnz(Qh)>0
                [i,j,k] = find(Qh);
                i = lmi_variables(i);
                j = lmi_variables(j);
                Q = sparse(i,j,k,nvars,nvars);
            end
        else
            % A relaxed problem should not calculate quadratic
            % decomposistion, fix!
            lmi_variables = getvariables(h);
            base = getbase(h);base= base(2:end);
            c=base(1)*zeros(nvars,1);
            c(lmi_variables) = base;
            Q = spalloc(nvars,nvars,0); 
            f = full(getbasematrix(h,0));
        end
    end
end
