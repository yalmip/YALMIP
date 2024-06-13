function p = detect_sdpmonotonicity(p)
p.sdpmonotinicity = [];
if p.options.cutsdp.nodefix && (p.K.s(1)>0)
    top = startofSDPCone(p.K);
    for i=1:length(p.K.s)
        n=p.K.s(i);
        for j=1:size(p.F_struc,2)-1
            X=full(reshape(p.F_struc(top:top+n^2-1,j+1),p.K.s(i),p.K.s(i)));
            X=(X+X')/2;
            v=real(eig(X+sqrt(eps)*eye(length(X))));
            if all(v>=0)
                sdpmonotinicity(i,j)=-1;
            elseif all(v<=0)
                sdpmonotinicity(i,j)=1;
            else
                sdpmonotinicity(i,j)=nan;
            end
        end
        top=top+n^2;
    end
else
    sdpmonotinicity=[];
end
