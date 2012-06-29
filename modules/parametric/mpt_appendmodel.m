function model = savemptmodel(model,Pfinal,Pn,Fi,Gi,details);

if length(Fi)>0
    if length(model) == 0
        model{1} = fakemptmodel(Pfinal, Pn, Fi, Gi, details.Ai, details.Bi, details.Ci);
        [H,K] = double(Pfinal);
        model{1}.epicost = [];
        model{1}.convex = 1;
    else
        anyqp = nnz([details.Ai{:}])>0;
        for i = 1:length(model)
            anyqp = anyqp | nnz([model{i}.Ai{:}])>0;
            if anyqp
                break
            end
        end
        if anyqp
            model = savemptmodelqp(model,Pfinal,Pn,Fi,Gi,details);
            return
        else
            newmodel = fakemptmodel(Pfinal, Pn, Fi, Gi, details.Ai, details.Bi, details.Ci);
            newmodel.epicost = [];
            newmodel.convex = 1;

            replace = zeros(length(model),1);
            discard = 0;
            for i = 1:length(model)

                Y = model{i}.Pfinal;
              %  quickeq( Pfinal, Y)
%                 if (Pfinal == Y) ~= quickeq( Pfinal, Y)
%                    1 
%                 end
                if  Pfinal == Y

                    if isempty(newmodel.epicost)
                        B = reshape([details.Bi{:}]',size(details.Bi{1},2),[])';
                        c = reshape([details.Ci{:}]',[],1);
                        [H,K] = double(Pfinal);
                        newmodel.epicost = generate_epicost(H,K,B,c);
                    end
                    if isempty(model{i}.epicost)
                        B = reshape([model{i}.Bi{:}]',size(model{i}.Bi{1},2),[])';
                        c = reshape([model{i}.Ci{:}]',[],1);
                        [H,K] = double(model{i}.Pfinal);
                        model{i}.epicost = generate_epicost(H,K,B,c);
                    end

                    if newmodel.epicost <= model{i}.epicost
                        discard = 1;
                        break
                    end

                    if newmodel.epicost >= model{i}.epicost
                        replace(i,1) = 1;
                    end

                elseif Pfinal >= model{i}.Pfinal

                    if isempty(newmodel.epicost)
                        B = reshape([details.Bi{:}]',size(details.Bi{1},2),[])';
                        c = reshape([details.Ci{:}]',[],1);
                        [H,K] = double(Pfinal);
                        newmodel.epicost = generate_epicost(H,K,B,c);
                    end
                    if isempty(model{i}.epicost)
                        B = reshape([model{i}.Bi{:}]',size(model{i}.Bi{1},2),[])';
                        c = reshape([model{i}.Ci{:}]',[],1);
                        [H,K] = double(model{i}.Pfinal);
                        model{i}.epicost = generate_epicost(H,K,B,c);
                    end
                    if newmodel.epicost >= model{i}.epicost
                        replace(i,1) = 1;
                    end

                elseif Pfinal <= model{i}.Pfinal

                    if isempty(newmodel.epicost)
                        B = reshape([details.Bi{:}]',size(details.Bi{1},2),[])';
                        c = reshape([details.Ci{:}]',[],1);
                        [H,K] = double(Pfinal);
                        newmodel.epicost = generate_epicost(H,K,B,c);
                    end
                    if isempty(model{i}.epicost)
                        B = reshape([model{i}.Bi{:}]',size(model{i}.Bi{1},2),[])';
                        c = reshape([model{i}.Ci{:}]',[],1);
                        [H,K] = double(model{i}.Pfinal);
                        model{i}.epicost = generate_epicost(H,K,B,c);
                    end

                    if newmodel.epicost <= model{i}.epicost
                        discard = 1;
                        break
                    end
                end
            end
            if ~discard
                model = {model{find(~replace)},newmodel};
            end
        end
    end
end

function model = savemptmodelqp(model,Pfinal,Pn,Fi,Gi,details);
newmodel = fakemptmodel(Pfinal, Pn, Fi, Gi, details.Ai, details.Bi, details.Ci);
newmodel.epicost = [];
newmodel.convex = 1;

replace = zeros(length(model),1);
discard = 0;

for i = 1:length(model)

    % Trivial pruning, why not...
    if quadraticLarger(details,model{i})
        if Pfinal == model{i}.Pfinal
            discard = 1;
            break;
        elseif Pfinal <= model{i}.Pfinal
            discard = 1;
            break;
        end
    end

    if Pfinal==model{i}.Pfinal

    elseif Pfinal >= model{i}.Pfinal

        doreplace = 1;
        for k = 1:length(details.Ai)
            Qnew = [details.Ai{k}   0.5*details.Bi{k}' ;0.5*details.Bi{k} details.Ci{k}];
            for j = 1:length(model{i}.Pfinal)
                Qold = [model{i}.Ai{j}  0.5*model{i}.Bi{j}';0.5*model{i}.Bi{j} model{i}.Ci{j}];
                if ~all(real(eig(full(Qold-Qnew)))>=-1e-12)
                    doreplace = 0;
                end
            end
        end
        if doreplace
            replace(i,1) = 1;
        end

    elseif Pfinal <= model{i}.Pfinal

        %                 if relaxationLarger(details,model{i},Pn,model{i}.Pn)
        %                 discard = 1;
        %                 break
        %             end
        %
        if length(Pn)==1
            discard = 1;
            Qnew = [details.Ai{1}   0.5*details.Bi{1}' ;0.5*details.Bi{1} details.Ci{1}];
            for j = 1:length(model{i}.Pfinal)
                Qold = [model{i}.Ai{j}  0.5*model{i}.Bi{j}';0.5*model{i}.Bi{j} model{i}.Ci{j}];
                if ~all(real(eig(full(Qnew-Qold)))>=-1e-12)
                    discard = 0;
                end
            end
            if discard
                break
            end
        end
%         
%         if length(Pn)==1
%             discard = 1;
%             Qnew = [details.Ai{1}   0.5*details.Bi{1}' ;0.5*details.Bi{1} details.Ci{1}];
%             for j = 1:length(model{i}.Pn)
%                 Qold = [model{i}.Ai{j}  0.5*model{i}.Bi{j}';0.5*model{i}.Bi{j} model{i}.Ci{j}];
%                 Pis = intersect(model{i}.Pn(j),newmodel.Pfinal);
%                 if isfulldim(Pis)
%                     [H,K] = double(Pis);
%                     if ~quadraticLarger2(Qnew,Qold,H,K)
%                         discard = 0;
%                     end
%                 end
%             end
%             if discard
%                 break
%             end
%         end
%         
    end
end

if ~discard
    model = {model{find(~replace)},newmodel};
end

function XbiggerY = relaxationLarger(X,Y,P1,P2)
XbiggerY = 1;
x = sdpvar(length(X.Ai{1}),1);
for k = 1:length(X.Ai)
    Xq = [X.Ai{k}   0.5*X.Bi{k}' ;0.5*X.Bi{k} X.Ci{k}];
    p1 = [x;1]'*Xq*[x;1];
    for j = 1:length(Y.Ai)
        Yq = [Y.Ai{j}   0.5*Y.Bi{j}' ;0.5*Y.Bi{j} Y.Ci{j}];
        if all(real(eig(full(Xq-Yq)))>0)
        else
            isc = intersect(P1(k),P2(j));
            if isfulldim(isc)
                p2 = [x;1]'*Yq*[x;1];
                [H,K] = double(isc);
                [xx,cc,vv]= solvemoment(set(H*x < K),p1-p2)
            end
        end
        %        if ~all(real(eig(Xq-Yq))>=-1e-12)
%            XbiggerY = 0;
%            return
%        end
    end
end


function XbiggerY = quadraticLarger2(Q1,Q2,H,K)

x = sdpvar(size(H,2),1);
obj = [x;1]'*Q1*[x;1]-[x;1]'*Q2*[x;1];

%sol = solvemoment(set(H*x <= K),obj,[],2);
sol = solvesdp(set(H*x <= K),obj,sdpsettings('solver','kktqp'));
%if relaxdouble(obj) > -1e-5
if double(obj) > -1e-5

    XbiggerY = 1;
else
    XbiggerY = 0;
end

function XbiggerY = quadraticLarger(X,Y)
XbiggerY = 1;
for k = 1:length(X.Ai)
    Xq = [X.Ai{k}   0.5*X.Bi{k}' ;0.5*X.Bi{k} X.Ci{k}];
    for j = 1:length(Y.Ai)
        Yq = [Y.Ai{j}   0.5*Y.Bi{j}' ;0.5*Y.Bi{j} Y.Ci{j}];
        if ~all(real(eig(full(Xq-Yq)))>=-1e-12)
            XbiggerY = 0;
            return
        end
    end
end

function epicost = generate_epicost(H,K,B,c)
epicost = polytope([B -ones(size(B,1),1);H zeros(size(H,1),1)],[-c;K]);



function C = fakemptmodel(Pfinal, Pn, Fi, Gi, Ai, Bi, Ci)

dummystruct.dummyfield = [];

C.sysStruct = dummystruct;
C.probStruct = dummystruct;
C.details.origSysStruct = dummystruct;
C.details.origProbStruct = dummystruct;
nr = length(Pn);

C.Pfinal = Pfinal;
C.Pn = Pn;
C.Fi = Fi;
C.Gi = Gi;
if isempty(Ai),
    Ai = cell(1, nr);
end
C.Ai = Ai;
C.Bi = Bi;
C.Ci = Ci;
C.dynamics = repmat(0, 1, nr);
C.overlaps = 0;



function status = quickeq(P,Q)

status = 1;
Q = struct(Q);
P = struct(P);
[ncP,nxP]=size(P.H);
[ncQ,nxQ]=size(Q.H);

if ncP ~=ncQ
    status = 0;
    return
end

Pbbox = P.bbox;
Qbbox = Q.bbox;
if ~isempty(Pbbox) & ~isempty(Qbbox),
    bbox_tol = 1e-4;
    if any(abs(Pbbox(:,1) - Qbbox(:,1)) > bbox_tol) | any(abs(Pbbox(:,2) - Qbbox(:,2)) > bbox_tol),
        % bounding boxes differ by more than abs_tol => polytopes cannot be equal
        status = 0;
        return
    end
    % we cannot reach any conclusion based solely on the fact that bounding
    % boxes are identical, therefore we continue...
end

status=1;
PAB=[P.H P.K];
QAB=[Q.H Q.K];

for ii=1:ncP
    %if all(sum(abs(QAB-repmat(PAB(ii,:),ncQ,1)),2)>abs_tol)
    %if all(sum(abs(QAB-PAB(ones(ncQ,1)*ii,:)),2)>abs_tol)
    Z = sum(abs(QAB-PAB(ones(ncQ,1)*ii,:)),2);
    if all(Z>1e-8)
        status=0;
        return;
    end
end

for ii=1:ncQ
    %if all(sum(abs(PAB-repmat(QAB(ii,:),ncP,1)),2)>abs_tol)
    %if all(sum(abs(PAB-QAB(ones(ncP,1)*ii,:)),2)>abs_tol)
    Z = sum(abs(PAB-QAB(ones(ncP,1)*ii,:)),2);
    if all(Z>1e-8)
        status=0;
        return;
    end
end
