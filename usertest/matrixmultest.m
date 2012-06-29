d1 = ceil(rand(1)*25+1);d2 = ceil(rand(1)*10+1);A = sprandn(d1,d2,0.1);Y =sdpvar(d2,d1);
tic;for i = 1:100;A*Y;end;toc