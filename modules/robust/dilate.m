function varargout = dilate(F,w,options)
% DILATE  Derives a matrix dilation
%
% [G,H,M] = DILATE(X,w,options) where X is a symmetric variable derives the
% decomposition and orthogonal complement used in a matrix dilation.
%
% X is decomposed as M(w)´G(x)M(w), and H(w) is an orthogonal complement to
% M(w), with affine dependence in w. These matrices can be used to apply
% the matrix dilation lemma to obtain a sufficient constraint affine in w
%
% M(w)´G(x)M(w) > 0  <== existence of W s.t with G + W*H' + H*W' > 0
%
% F = DILATE(F,w,options) where F is a SET is used to construct the dilated
% uncertain constraint, i.e an SDP constraint where polynomial dependence
% w.r.t uncertain variables in all SDP constraints in F are converted
% (conservatively) to affine dependence using the matrix dilation approach. 
%
% The SDPSETTING options structure is optional, but recommended if you call
% dilate repeatedly (avoids calling SDPSETTINGS inside DILATE)
%
% See Yasuaki OISHI, A Region-Dividing Approach to Robust Semidefinite
% Programming and Its Error Bound,DEPARTMENT OF MATHEMATICAL INFORMATICS
% GRADUATE SCHOOL OF INFORMATION SCIENCE AND TECHNOLOGY THE UNIVERSITY OF
% TOKYO BUNKYO-KU, TOKYO 113-8656, JAPAN , February 2006
%
% J. Löfberg, Improved matrix dilations for robust semidefinite programming
% Division of Automatic control, Department of Electrical Engineering, 
% Technical report 2753, Linköpings universitet, November 2006
%
% See also ROBUSTIFY, SOLVEROBUST, UNCERTAIN

% Author Johan Löfberg
% $Id: dilate.m,v 1.9 2009-04-29 11:38:00 joloef Exp $

if nargin < 3
    options = sdpsettings;
end

if isa(F,'constraint')
    F = set(F);
end

if isa(F,'sdpvar')
    [G,H,M] = matrix_dilate(F,w,options);
    varargout{1} = G;
    varargout{2} = H;
    varargout{3} = M;
elseif isa(F,'lmi')
    Fnew = [];
    if nargin == 1
        w = [];
    else
        w = w(:);
    end
    unc_declarations = is(F,'uncertain');
    if any(unc_declarations)
        w = [w;recover(getvariables(sdpvar(F(find(unc_declarations)))))];
    end
    for i = 1:length(F)
        if (is(F(i),'sdp') | ((length(sdpvar(F(i))) == 1) & is(F(i),'elementwise'))) & max(degree(sdpvar(F(i)),w))>0
            [G,H,M] = matrix_dilate(sdpvar(F(i)),w,options);
            W = sdpvar(size(G,1),size(H,2));
            Fnew = Fnew + set(G + W*H' + H*W' >= 0);
        else
            Fnew = Fnew + F(i);
        end
    end
    varargout{1} = Fnew;
end


function [G,H,M] = matrix_dilate(F,w,options)
% Given an SDP F(x,w) with F polynomial in w, DILATE rewrites the problem
% to F(x,w)=M(w)'G(x)M(w), and derives the orhogonal complements H(w) to
% M(w), to be used in the dilated constraint G + W*H' + H*W'

% Parametric variables in SOS i.e. the certain variables
x = recover(setdiff(depends(F),depends(w)));
v = monolist(w,ceil(max(degree(F,w)/2)));
options.sos.newton = 0;
options.sos.model = 2;
options.sos.scale = 0;
options.verbose = 0;
options.sos.conggruence = 0;
[G,dummy,m] = compilesos(set(sos(F)),[],options,x,v);
G = sdpvar(G);

allvariables = getvariables(F);
wvariables   = getvariables(w);
Fbasis = getbase(F);

% Make sure it is ordered according to internal index
w = recover(wvariables);

% Degrees w.r.t the uncertain variables
d = degree(recover(allvariables),w);

% robustify code by removing unused variables
w = w(find(d));
wvariables =  wvariables(find(d));

% Degrees w.r.t the uncertain variables
monomtable = yalmip('monomtable');
d = max(sum(monomtable(allvariables,wvariables),2));
n = size(F,1);

% Sufficiently many monomials
v = m{1};%monolist(w,max(d));
v = v(1:(size(v,1)/size(v,2)),1);
% Some numerical format on these variables
for i = 2:length(v)
    vvariables(i,1) = getvariables(v(i));  
end
monomtable = yalmip('monomtable');
wmonoms = [zeros(1,length(w));monomtable(vvariables(2:end),wvariables)];

% Find the linear certain term in the matrix
linear_indicies = find(sum(monomtable(allvariables,wvariables),2) ==0);
if length(linear_indicies) > 0
    F0 = reshape(Fbasis(:,1) + Fbasis(:,1+linear_indicies)*recover(allvariables(linear_indicies)),n,n);
else
    F0 = reshape(Fbasis(:,1),n,n);
end

% now find the matrix that multiplies with each monomial in M

% Fi = [];
% Fmonoms = monomtable(allvariables,wvariables);
% for i = 2:length(v)    
%     vmonoms = monomtable(vvariables(i),wvariables);
%     index = findrows(Fmonoms,vmonoms);
%     if isempty(index)
%         Fi = [Fi zeros(n)];
%     else
%         temp = monomtable(allvariables(index),:);
%         temp(wvariables) = 0;
%         base = reshape(Fbasis(:,1+index),n,n);
%         if nnz(temp) == 0
%             Fi = [Fi base];
%         else
%             Fi =  [Fi base*recover(allvariables(find(temp)))];
%         end        
%     end
% end
% G = [F0 Fi/2;Fi'/2 zeros(n*(length(v)-1))];

% The outer factor
M = m{1};%kron(v,eye(n));

% Now create an orthogonal complement to M
ii = [];
jj = [];
ss = [];
for i = 2:length(v)
    monom = wmonoms(i,:);
    [dummy,index] = max(monom);
    monomnew = monom;
    monomnew(index) = monomnew(index) - 1;   
    ii = [ii i];
    jj = [jj i-1];
    ss = [ss 1];
    ii = [ii findrows(wmonoms,monomnew)];
    jj = [jj i-1];
    ss = [ss -recover(wvariables(index))];
end
H = kron(sparse(ii,jj,ss),eye(n));
