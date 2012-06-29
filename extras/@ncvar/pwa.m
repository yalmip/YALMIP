function [p,Bi,Ci,Pn,Pfinal] = PWA(h,Xdomain)
% PWA Tries to create a PWA description
%
% [p,Bi,Ci,Pn,Pfinal] = PWA(h,X)
%
% Input
%  h  : scalar SDPVAR object
%  X  : SET object
%
% Output
%
%  p  : scalar SDPVAR object representing the PWA function
%  Bi,Ci,Pn,Pfinal : Data in MPT format
%
% The command tries to expand the nonlinear operators
% (min,max,abs,...) used in the variable h, in order
% to generate an epi-graph model. Given this epigraph model,
% it is projected to the variables of interest and the
% defining facets of the PWA function is extracted. The second 
% argument can be used to limit the domain of the PWA function.
% If no second argument is supplied, the PWA function is created
% over the domain -100 to 100.
% 
% A new sdpvar object p is created, representing the same 
% function as h, but in a slightly different internal format.
% Additionally, the PWA description in MPT format is created.
%
% The function is mainly inteded to be used for easy 
% plotting of convex PWA functions

% Author Johan Löfberg 
% $Id: pwa.m,v 1.1 2006-08-10 18:00:22 joloef Exp $   


t = sdpvar(1,1);
[F,failure,cause] = expandmodel(set(h<t),[],sdpsettings('allowmilp',0));
if failure
    error(['Could not expand model (' cause ')']);
    return
end

% Figure out what the actual original variables are
% note, by construction, they 
all_initial = getvariables(h);
all_extended = yalmip('extvariables');
all_variables = getvariables(F);
gen_here = getvariables(t);
non_ext_in = setdiff(all_initial,all_extended);
lifted = all_variables(all_variables>gen_here);
vars = union(setdiff(setdiff(setdiff(all_variables,all_extended),gen_here),lifted),non_ext_in);
nx = length(vars);

X = recover(vars);

if nargin == 1
    Xdomain = set(-100 < X < 100);
else
    Xdomain = set(-10000 < X < 10000)+Xdomain;
end

[Ai,Bi,Ci,Pn] = generate_pwa(F,t,X,Xdomain,nx);

Pfinal = union(Pn);
sol.Pn = Pn;
sol.Bi = Bi;
sol.Ci = Ci;
sol.Ai = Ai;
sol.Pfinal = Pfinal;
p = pwf(sol,X,'convex'); 
% 
% binarys = recover(all_variables(find(ismember(all_variables,yalmip('binvariables')))))
% if length(binarys) > 0
%     
%     Binary_Equalities = [];
%     Binary_Inequalities = [];
%     Mixed_Equalities = [];
%     top = 1;
%     for i = 1:length(F)
%         Fi = sdpvar(F(i));
%         if is(F(i),'equality')
%             if all(ismember(getvariables(Fi),yalmip('binvariables')))
%                 Binary_Equalities = [Binary_Equalities;(top:top-1+prod(size(Fi)))'];
%                 Mixed_Equalities = [Mixed_Equalities;(top:top-1+prod(size(Fi)))'];
%             end
%         else
%             if all(ismember(getvariables(Fi),yalmip('binvariables')))
%                 Binary_Inequalities = [Binary_Inequalities;(top:top-1+prod(size(Fi)))'];
%             end
%         end
%         top = top+prod(size(Fi))';
%     end
%     P = sdpvar(F);
%     P_ineq = extsubsref(P,setdiff(1:length(P),[Binary_Equalities; Binary_Inequalities]))
%     P_binary_eq = extsubsref(P,Binary_Equalities);HK1 = getbase(P_binary_eq);
%     P_binary_ineq = extsubsref(P,Binary_Inequalities);HK2 = getbase(P_binary_ineq);
%     nbin = length(binarys);
%     enums = dec2decbin(0:2^nbin-1,nbin)'
%     if isempty(HK2)
%         HK2 = HK1*0;
%     end
%     for i = 1:size(enums,2)
%         if all(HK1*[1;enums(:,i)]==0)
%             if all(HK2*[1;enums(:,i)]>=0)
%                 Pi = replace(P_ineq,binarys,enums(:,i))
%             end
%         end
%     end    
%     
% else
%     [Ai,Bi,Ci,Pn] = generate_pwa(F,t,X,Xdomain,nx);
% end
% 



function [Ai,Bi,Ci,Pn] = generate_pwa(F,t,X,Xdomain,nx)

% Project, but remember that we already expanded the constraints
P = polytope(projection(F+set(t<10000)+Xdomain,[X;t],[],1));Xdomain = polytope(Xdomain);
[H,K] = double(P);
facets = find(H(:,end)<0);

region = find(~H(:,end) & any(H(:,1:nx),2) );
Hr = H(region,1:nx);
Kr = K(region,:);

H = H(facets,:);
K = K(facets);
K = K./H(:,end);
H = H./repmat(H(:,end),1,size(H,2));

nx = length(X);
Pn = [];
cib = [H(:,1:nx) K];
Ai = {};
Bi = cell(0);
Ci = cell(0);
if length(Kr > 0)
    Xdomain = intersect(Xdomain,polytope(Hr,Kr));
end
[Hr,Kr] = double(Xdomain);

for i = 1:length(K)
    j = setdiff(1:length(K),i);
    HiKi = repmat(cib(i,:),length(K)-1,1)-cib(j,:);
    Pi = polytope([HiKi(:,1:nx);Hr],[HiKi(:,end);Kr]);

    if isfulldim(Pi)
        Pn = [Pn Pi];
        Bi{end+1} = -cib(i,1:end-1);
        Ci{end+1} = cib(i,end);
        Ai{end+1} = [];
    end
end