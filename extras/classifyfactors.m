function [constants,general,singles,pairs] = classifyfactors(L,M,R)

nfac = length(L);
used = zeros(1,nfac);

constants = 0;
general = 0;
singles = {};
pairs = {};
for i = 1:length(L)
    if isa(M{i},'double')
        constants = constants + L{i}*M{i}*R{i};
        used(i)=1;
    end
end

for i = 1:length(L)
    if ~used(i)      
        temp = struct(M{i});
        if ~(is(M{i},'sdpcone') || is(M{i},'lpcone') || isequal(temp.originalbasis,'skew') || isequal(temp.originalbasis,'diagonal'))
            general = general + L{i}*M{i}*R{i};
            used(i) = 1;
        elseif size(M{i},2)~=size(R{i},1) | size(M{i},1)~=size(L{i},2)
            % Something like scalar*matrix
            % simplifies code if we treat that as a general
            general = general + L{i}*M{i}*R{i};
            used(i) = 1;
        end
    end
end

% sameL = [];
% sameR = [];
% for i = 1:length(L)-1
%     for j = i+1:length(L)
%         if isequal(M{i},M{j})
%             if isequal(L{i},L{j})
%                 mmm(i,j) = 1;
%                 mmm(j,i) = 1;
%                 sameL = [sameL;i j];
%             elseif isequal(R{i},R{j})
%                 sameR = [sameR;i j];
%                 mmm(i,j) = 1;
%                 mmm(j,i) = 1;
%             end            
%         end
%     end
% end
% if length(sameL)>1
%     for i = 1:size(sameL,1)-1
%         for j = 2:size(sameL,2)
%             AA = sameL([i j],:);
%             BB = sameR([i j],:);
%             if AA == BB'
%                 candidate = unique([AA(:);BB(:)]);
%             end
%         end
%     end
% end

for i = 1:length(L)
    
     if issymmetric(M{i}) && isequal(find(L{i}),find(R{i}'))
        Rttmp = R{i}';
        tmp = L{i}(L{i}~=0)./(Rttmp(Rttmp~=0));
     else
         tmp = [];
     end
    
    if ~used(i)
        if issymmetric(M{i}) && (isequal(L{i},R{i}')) 
            singles{end+1}.L = L{i};
            singles{end}.M = M{i};
            singles{end}.R = R{i};
            singles{end}.negated = 0;
            used(i)=1;
        elseif issymmetric(M{i}) && (isequal(L{i},-R{i}')) 
            singles{end+1}.L = L{i};
            singles{end}.M = M{i};
            singles{end}.R = R{i};
            singles{end}.negated = 1;
            used(i)=1;          
        elseif issymmetric(M{i}) && isequal(find(L{i}),find(R{i}')) &  length(unique(tmp)) == 1       
        %    Rttmp = R{i}';
        %    tmp = L{i}(L{i}~=0)./(Rttmp(Rttmp~=0));
            
         %   if length(unique(tmp)) == 1                   
                  singles{end+1}.L = L{i}/sqrt(abs(tmp(1)));
                  singles{end}.M = M{i};
                  singles{end}.R = R{i}*sqrt(abs(tmp(1)));
                  if tmp(1) > 0
                      singles{end}.negated = 0;
                  else
                      singles{end}.negated = 1;
                      
                  end
                  used(i)=1;                
        %  end
    
        else            
            for j = i+1:length(L)
                if isa(L{i}*M{i}*R{i}-R{j}'*M{j}'*L{j}','double')
                    pairs{end+1}.L = L{i};
                    pairs{end}.M = M{i};
                    pairs{end}.R = R{i};
                    used(j)=1;
                    used(i)=1;
                    break
                end
            end
        end
    end
end


for k = find(~used)
    if isequal(size(M{k}),[1 1])
        general = general + L{k}*M{k}*R{k};
    else
        error('Couldn''t classify all factors');
    end
end


