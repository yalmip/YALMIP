function   model = yalmip2pensdp(interfacedata);

penstruct = sedumi2pen(interfacedata.F_struc,interfacedata.K,interfacedata.c,interfacedata.x0);
ops = struct2cell(interfacedata.options.pensdp);ops = [ops{1:end}];
penstruct.ioptions = ops(1:8);
penstruct.foptions = ops(9:end);
penstruct.ioptions(4) = interfacedata.options.verbose;
penstruct.ioptions = penstruct.ioptions;
penstruct.foptions = penstruct.foptions;
if penstruct.mconstr == 0
    penstruct.msizes = [];
end
model.penstruct = penstruct;