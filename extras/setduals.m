function setduals(F,D_struc,K)

D_struc = full(D_struc);
if ~isempty(D_struc)
   
    % Equality constraints might have been converted
    % to inequalities
    if isfield(K,'fold')
       D_struc = [-D_struc(1:K.fold)+D_struc(1+K.fold:2*K.fold);D_struc(2*K.fold+1:end)];
       K.f = K.fold;
       K.l = K.l-2*K.fold;
    end
    
    Z = dual2cell(D_struc,K);
    
    Finfo = lmiinfo(F);
    lmi_index = [];
    j=1;
    
    if ~isempty(Finfo.sdp)
        for i = 1:size(Finfo.sdp,1)
            lmi_index = [lmi_index;Finfo.sdp(i,3)];
            duals{j}=(Z.s{i}+Z.s{i}')/2;j = j+1;
        end
    end
    
    if ~isempty(Finfo.soc)
        for i = 1:size(Finfo.soc,1)
            lmi_index = [lmi_index;Finfo.soc(i,3)];
            duals{j}=Z.q{i};j = j+1;
        end
    end
    
    if ~isempty(Finfo.rlc)
        for i = 1:size(Finfo.rlc,1)
            lmi_index = [lmi_index;Finfo.rlc(i,3)];
            duals{j}=Z.r{i};j = j+1;
        end
    end
    
    if ~isempty(Finfo.lin)
        top=1;
        for i = 1:size(Finfo.lin,1)
            lmi_index = [lmi_index;Finfo.lin(i,3)];
            duals{j}=reshape(Z.l(top:top+Finfo.lin(i,1)*Finfo.lin(i,2)-1),Finfo.lin(i,1),Finfo.lin(i,2));j = j+1;
            top = top+Finfo.lin(i,1)*Finfo.lin(i,2);
        end
    end
    if ~isempty(Finfo.equ)
        top=1;
        for i = 1:size(Finfo.equ,1)
            lmi_index = [lmi_index;Finfo.equ(i,3)];
            duals{j}=reshape(Z.f(top:top+Finfo.equ(i,1)*Finfo.equ(i,2)-1),Finfo.equ(i,1),Finfo.equ(i,2));j = j+1;
            top = top+Finfo.equ(i,1)*Finfo.equ(i,2);
        end
    end
    
    if ~isempty(lmi_index)
        yalmip('setdual',lmi_index,duals);
    end
    
end
