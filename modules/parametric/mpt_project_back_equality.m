function [Fi,Gi,details] = mpt_project_back_equality(Matrices,Fi,Gi,details,OriginalMatrices);
if isempty(Matrices.getback)
    return
else
    %Original solution given by S1*z+S2*x+S3
    S1 = Matrices.getback.S1;
    S2 = Matrices.getback.S2;
    S3 = Matrices.getback.S3;
    for i=1:length(Fi)
        if isempty(Fi{i})
            Fi{i} = S2;
            Gi{i} = S3;
        else
        Fi{i} = S1*Fi{i} + S2;
        Gi{i} = S1*Gi{i} + S3;
        end
    end
    if nargin>2
        if OriginalMatrices.qp
            for i=1:length(Fi)
                % FIX : Check this...
                details.Ai{i} = 0.5*Fi{i}'*OriginalMatrices.H*Fi{i} + 0.5*(OriginalMatrices.F*Fi{i}+Fi{i}'*OriginalMatrices.F') + OriginalMatrices.Y;
                details.Bi{i} = OriginalMatrices.Cf*Fi{i}+Gi{i}'*OriginalMatrices.F' + Gi{i}'*OriginalMatrices.H*Fi{i} + OriginalMatrices.Cx;
                details.Ci{i} = OriginalMatrices.Cf*Gi{i}+0.5*Gi{i}'*OriginalMatrices.H*Gi{i} + OriginalMatrices.Cc;
            end
        else
            for i=1:length(Fi)
                details.Bi{i} = OriginalMatrices.H*Fi{i} + OriginalMatrices.F;
                details.Ci{i} = OriginalMatrices.H*Gi{i};
            end
        end
    end
end
