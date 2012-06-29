function [BlockedQ,residuals] = postprocesssos(BlockedA,Blockedb,BlockedQ,sparsityPattern,options);

BlockedQ=applysparsity(BlockedQ,sparsityPattern);

for passes = 1:1:options.sos.postprocess
    for constraint = 1:length(BlockedQ)
        mismatch = computeresiduals(BlockedA,Blockedb,BlockedQ,constraint);
        [ii,jj ]= sort(abs(mismatch));
        jj = flipud(jj);
        for j = jj(:)'%1:size(BlockedA{constraint}{1},1)
            if abs(mismatch(j))>0
                for i = 1:length(BlockedA{constraint})
                    n=sqrt(size(BlockedA{constraint}{i},2));
                    Ai=reshape(BlockedA{constraint}{i}(j,:),n,n);
                    nnzAi = nnz(Ai);
                    if nnzAi>0
                        dAi = Ai*mismatch(j)/nnzAi;
                        Qi = BlockedQ{constraint}{i};
                        %                        [R,p] = chol(BlockedQ{constraint}{i}-Ai*mismatch(j)/nnzAi);
                        [R,p] = chol(Qi-dAi);
                        if p
                            %                           gevps=eig(BlockedQ{constraint}{i},full(Ai)*mismatch(j)/nnzAi);
                            gevps=eig(Qi,full(dAi));
                            gevps=gevps(gevps>=0);
                            gevps=gevps(~isinf(gevps));
                            if isempty(gevps)
                                gevps=1;
                            end
                            lambda=max(0,min(1,min(gevps)));
                            %                            [R,p] = chol(BlockedQ{constraint}{i}-Ai*lambda*mismatch(j)/nnzAi);
                            [R,p] = chol(Qi-dAi*lambda);
                        else
                            lambda = 1;
                        end
                        %                        dAi = Ai*mismatch(j)/nnzAi;
                        if ~p
                            %                            BlockedQ{constraint}{i}=BlockedQ{constraint}{i}-Ai*lambda*mismatch(j)/nnzAi;
                            BlockedQ{constraint}{i}=Qi-dAi*lambda;
                            mismatch(j)=mismatch(j)*(1-lambda);
                        else
                            %lambda=1;
                            while lambda>1e-4 & (p~=0)
                                lambda=lambda/sqrt(2);
                                %                                [R,p]= chol(BlockedQ{constraint}{i}-Ai*lambda*mismatch(j)/nnzAi);
                                [R,p]= chol(Qi-dAi*lambda);
                            end
                            %                            if min(eig(BlockedQ{constraint}{i}-Ai*lambda*mismatch(j)/nnzAi))>=0
                            if min(eig(Qi-dAi*lambda))>=0
                                %    BlockedQ{constraint}{i}=BlockedQ{constraint}{i}-Ai*lambda*mismatch(j)/nnzAi;
                                BlockedQ{constraint}{i}=Qi-dAi*lambda;
                                mismatch(j)=mismatch(j)*(1-lambda);
                            end
                        end
                    end
                end
            end
        end
    end
end

for constraint = 1:length(BlockedQ)
    residuals(constraint,1) = norm(computeresiduals(BlockedA,Blockedb,BlockedQ,constraint),'inf');
end

function  mismatch = computeresiduals(BlockedA,Blockedb,BlockedQ,constraint);
lhs=0;
for k=1:length(BlockedA{constraint})
    lhs=lhs+BlockedA{constraint}{k}*BlockedQ{constraint}{k}(:);
end
mismatch = lhs-Blockedb{constraint};

function BlockedQ=applysparsity(BlockedQ,sparsityPattern);
if ~isempty(sparsityPattern)

    for i = 1:length(BlockedQ)
        for j =  1:length(BlockedQ{i})
            BlockedQ{i}{j}(sparsityPattern{i}{j}) = 0;
        end
    end
end

