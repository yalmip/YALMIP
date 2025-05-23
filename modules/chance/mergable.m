function yesno = mergable(a,b)
%normalvariants = {'normal','normalf','normalm','mvnrnd','mvrndfactor'};
normalvariants = {'normal','mvnrnd','mvrndfactor'};
if isequal(a,b) && ~isequal(a,'dro')
    yesno = 1;   
    % normal, scalar elementwise normal (parameterized in std. dev.)
    % normalf, obsolete
    % normalm, obsolete
    % mvrnd, multivariable normal, (parameterized in covariance)
    % mvrndfactor, multivariable normal, (parameterized in factor covariance)
elseif any(strcmp(a,normalvariants)) && any(strcmp(b,normalvariants))
    yesno = 1;
else
    yesno = 0;
end