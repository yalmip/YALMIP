function new_x = monolist(x,dmax,dmin)
% MONOLIST Generate monomials
%
% y = MONOLIST(x,dmax,dmin)
%
% Returns the monomials [1 x(1) x(1)^2 ... x(1)^dmax(1) x(2) x(1)x(2)  etc...]
%
% >>sdpvar x y z
% >>sdisplay(monolist([x y z],4))
%
%  Input
%     x      : Vector with SDPVAR variables
%     dmax   : Integers > 0
%
%  See also POLYNOMIAL, DEGREE

% Flatten
x = reshape(x,1,length(x));
x_orig = x;

if nargin == 3
    if length(dmin)>1 || any(dmin > dmax) || ~isequal(dmin,fix(dmin))
        error('dmin has to be an integer scalar larger than dmax');
    end
elseif nargin == 2
    dmin = 0;
end

if (length(dmax)==1 | all(dmax(1)==dmax)) & islinear(x) & ~isa(x,'ncvar')
    dmax = dmax(1);
    % Powers in monomials
    powers = monpowers(length(x),dmax);
   
    powers = powers(find(sum(powers,2)>=dmin),:);
    
    % Use fast method for x^alpha
    if isequal(getbase(x),[zeros(length(x),1) eye(length(x))])
        new_x = recovermonoms(powers,x);
        return
    end

    % Monolist on dense vectors is currently extremely slow, but also
    % needed in some applications (stability analysis using SOS) For
    % performance issue, the code below is hard-coded for special cases
    % FIX : Urgent, find underlying indexing...
    
    % Vectorize quadratic and quadrtic case
    if dmax==2 & length(x)>1
        V=x.'*[1 x];
        ind=funkyindicies(length(x));
        new_x = [1 V(ind(:)).'].';
        return
    elseif (length(x)==4 & dmax==6)
        
         ind =[    1           2           3           4           5           6           7           8,

           9          10          11          12          13          14          15          16,

          17          18          19          20          21          22          23          24,

          25          26          27          28          29          30          31          32,

          33          34          49          50          51          52          86          53,

          54          55          89          56          57          91          58          92,

         126          59          60          61          95          62          63          97,

          64          98         132          65          66         100          67         101,


         135          68         102         136         170         185         186         187,

         188         222         256         189         190         191         225         259,


         192         193         227         261         194         228         262         296,

         330         364         195         196         197         231         265         198,

         199         233         267         200         234         268         302         336,

         370         201         202         236         270         203         237         271,

         305         339         373         204         238         272         306         340,

         374         408         442         476         510         525         526         527,

         528         562         596         630         529         530         531         565,


         599         633         532         533         567         601         635         534,


         568         602         636         670         704         738         772         806,


         840         535         536         537         571         605         639         538,

         539         573         607         641         540         574         608         642,

         676         710         744         778         812         846         541         542,

         576         610         644         543         577         611         645         679,

         713         747         781         815         849         544         578         612,


         646         680         714         748         782         816         850         884,


         918         952         986        1020        1054        1088        1122        1156,

        1190          0            0          0          0             0           0            0];
        ind = ind';
        ind = ind(find(ind));
        v=monolist(x,3);
        V=v(2:end)*v.';
        new_x = [1;V(ind(:))];
        return
    elseif dmax==4 & (1<=length(x)) & length(x)<=4 %& length(x)>1
        v=monolist(x,2);
        V=v(2:end)*v.';

        % Cone to generate indicies
        %p = sdpvar(n,1);
        %v = monolist(p,2);
        %V = v(2:end)*v';V=V(:);
        %m = monolist(p,4)
        %ind = [];
        %for i = 2:length(m)
        % ind = [ind min(find(~any(V-m(i))))];
        %end

        switch length(x)
            case 1
                new_x = [1; V([1 2 4 6]')];
                return
            case 2
                new_x = [1;V([1 2 3 4 5 8 9 10 15 18 19 20 25 30]')];
                return;
            case 3
                new_x=[1;V([1 2 3 4 5 6 7 8 9 13 14 15 24 16 17 26 18 27 36 40 41 42 51 60 43 44 53 62 45 54 63 72 81 90]')];
                return
            case 4
                new_x=[1;V([    1     2     3     4     5     6     7     8     9    10    11    12    13    14    19    20    21    35    22    23    37    24    38    52    25    26    40    27    41    55    28    42 56    70    75    76    77    91   105    78    79    93   107    80    94   108   122   136  150    81    82    96   110    83    97   111   125   139   153    84    98   112   126   140  154   168   182   196   210]')];
                return
            otherwise
        end
    end


    % Na, we have things like (c'x)^alpha
    % precalc x^p
    for i = 1:length(x)
        temp = x(i);
        precalc{i,1} = temp;
        for j = 2:1:dmax
            temp = temp*x(i);
            precalc{i,j} = temp;
        end
    end

    new_x = [];

    for i = 1:size(powers,1) % All monomials
        temp = 1;

        for j = 1:size(powers,2) % All variables
            if powers(i,j)>0
                temp = temp*precalc{j,powers(i,j)};
            end
        end
        new_x = [new_x temp];

    end

else

    dmax = dmax(:)*ones(1,length(x_orig));

    x = [1 x];

    % Lame loop to generate all combinations
    new_x = 1;
    for j = 1:1:max(dmax)
        temp = [];
        for i = 1:length(x)
            temp = [temp x(i)*new_x];
        end
        new_x = temp;
        new_x = fliplr(unique(new_x));
        new_degrees = degree(new_x,x(2:end));
        remv = [];
        for i = 1:length(dmax);
            if new_degrees(i)>=dmax(i)
                x = recover(setdiff(getvariables(x),getvariables(x_orig(i))));
                x = [1;x(:)];
                remv = [remv i];
            end
        end
        dmax = dmax(setdiff(1:length(dmax),remv));
    end
end
new_x = reshape(new_x(:),length(new_x),1);
if dmin > 0
    for i = 1:length(new_x)
        if sum(powers(i,:)) < dmin
            keep(i) = 0;
        else
            keep(i) = 1;
        end
    end
    if any(keep==0)
        new_x = new_x(find(keep));
    end
end

function ind = funkyindicies(n)

M=reshape(1:n*(n+1),n,n+1);
ind = M(:,1)';
for i = 1:n
    ind = [ind M(i,2:i+1)];
end


