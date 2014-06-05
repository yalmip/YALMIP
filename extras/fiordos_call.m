function  x = fiordos_call(solver,param,B0,b0,mask,map,dimout)

mparams.be = b0 + B0*param(mask); 
mres = solver(mparams);
allx = [nan;mres.x];% map=0 corresponds to stuff that really weren't in the model
x = reshape(allx(1+map),dimout);

