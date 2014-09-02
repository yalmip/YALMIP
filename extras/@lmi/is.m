function YESNO = is(F,property)
%IS   Check property of constraint.
%   d = IS(x,property) returns 1 if 'property' holds
%
%   Properties possible to test are: 'elementwise', 'sdp', 
%   'socc', 'equality', 'lmi', 'linear', 'kyp', 'sos'

if length(F.LMIid)==0
    YESNO = 0;
else
    
    F = flatten(F);
      
    YESNO=zeros(length(F.LMIid),1);
    switch property
        case 'dualized'
            YESNO = F.dualized == 1;
            
        case 'chance'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = ~isempty(Fi.confidencelevel);
            end
            
           case 'meta'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==56);
            end
        case 'equality'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==3);
            end
        case {'element-wise','elementwise'}
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==2);
            end
        case {'socc','socp'}
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==4);
            end
        case {'vecsocp'}
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==54);
            end            
        case 'pcone'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==20);
            end    
        
        case 'sdpcone'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = is(Fi.data,'sdpcone');
            end
            
        case 'realsdpcone'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = is(Fi.data,'realsdpcone');
            end
        
        case 'complexsdpcone'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = is(Fi.data,'complexsdpcone');
            end                
        case 'sdp'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = Fi.type==1;
            end    
        case 'lmi'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = full((((Fi.type==1) | (Fi.type==9)) | ((Fi.type==2) & (prod(size(Fi.data))==1))) & (islinear(Fi.data)));
            end
        case 'linear'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = islinear(Fi.data);
            end
        case 'kyp'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==9);
            end
        case 'sos2'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==50);
            end
        case 'sos1'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==51);
            end            
        case 'semivar'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==52);
            end            
        case 'semiintvar'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==53);
            end               
       case 'complementarity'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==55);
            end     
            
        case 'sos'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = (Fi.type==11);
            end
        case 'eig'
            YESNO(i,1) = (Fi.type==10);
        case 'sigmonial'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                monomtable = yalmip('monomtable');
                monomtable = monomtable(getvariables(Fi.data),:);
                YESNO(i,1) = any(find(any(0>monomtable,2) | any(monomtable-fix(monomtable),2)));
            end
        case 'binary'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = Fi.type ==  8;
            end
        case 'integer'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = Fi.type ==  7;
            end
       case 'parametric'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = Fi.type ==  13;
            end            
        case 'uncertain'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = Fi.type ==  15;
            end
        case 'random'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = Fi.type ==  16;
            end
            
            
        case 'logic'
            allextvars = yalmip('extvariables');
            if isempty(allextvars)
                YESNO(i,1) = 0;
            else
                for i = 1:length(F.clauses)
                    Fi = F.clauses{i};
                    xi = getvariables(Fi.data);
                    lgc = find(ismembc(xi,allextvars));
                    if ~isempty(lgc)
                        for j = lgc
                            variable = xi(j);
                            extstruct = yalmip('extstruct',11);
                            if isequal(extstruct.fcn,'or') | isequal(extstruct.fcn,'and')
                                YESNO(i,1) = 1;
                                break
                            end
                        end
                    end
                end                
            end
        case 'lowrank'            
             for i = 1:length(F.clauses)                
                  Fi = F.clauses{i};
                YESNO(i,1) = Fi.type == 14;
               % YESNO(i,1) = is(Fi.data,'complex');
             end                        
        case 'complex'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = isa(Fi.data,'sdpvar') && is(Fi.data,'complex');
            end
        case 'interval'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = isa(Fi.data,'sdpvar') && is(Fi.data,'interval');
            end            
        case 'real'
            for i = 1:length(F.clauses)
                Fi = F.clauses{i};
                YESNO(i,1) = isa(Fi.data,'sdpvar') && is(Fi.data,'real');
            end
        otherwise
            YESNO = error('Huh?');
    end
    %  end
end

YESNO = full(YESNO);