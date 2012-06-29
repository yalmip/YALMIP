function setsos(X,h,ParametricVariables,Q,v)
%SETSOS Internal function

yalmip('setsos',X.extra.sosid,h,ParametricVariables,Q,v);
