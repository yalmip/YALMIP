function model = droChanceFilter(b,c,distribution,gamma,w,options)
% Chance filter for distribution only specified by mean and variance
%themean    = distribution.parameters{2};
%covariance = distribution.parameters{3};
%alpha = distribution.parameters{4};
alpha = .005;
themean = [0 0];
covariance = 1/3;
h =  -c;
h0 = -b;

beta = sdpvar(1);
tau = sdpvar(1);
s = sdpvar(1);
z = sdpvar(1);
v = sdpvar(1);
n1 = sdpvar(1);
n2 = sdpvar(1);
if isinf(alpha)
    model = [z >= 0,
        v >= 0, n1 >= 0, n2 >= 0,
        beta+(1/gamma)*(v+s) <= 0,
        z+s >= norm([z-s;tau;h*covariance]),
        h == n1-n2,
        v-h0+beta+tau-(n1*themean(2)-n2*themean(1))-z >= 0,
        ];
else
    model = [z >= 0,
        v >= 0, n1 >= 0, n2 >= 0,
        beta+(1/gamma)*(v+s) <= 0,
        z+s >= norm([z-s;tau;sqrt((alpha+2)/alpha)*h*covariance]),
        h == n1-n2,
        v-h0+beta+(alpha/(alpha+1))*tau-(n1*themean(2)-n2*themean(1))-(alpha/(alpha+1))^2*z >= 0,
        ];
end