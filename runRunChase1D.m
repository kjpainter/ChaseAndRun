params.domainlength    = 10;
params.gridCells       = 20;

params.BCs = 'pp'; % periodic BCs

params.alpha_n1n1 = 3;
params.alpha_n1n2 =-3;
params.alpha_n2n1 = 3;
params.alpha_n2n2 = 3;


params.xi_n1n1 = 1.0;
params.xi_n1n2 = 1.0;
params.xi_n2n1 = 1.5;
params.xi_n2n2 = 1.5;

params.r = 0.0;
params.delta = 1.0;

params.fswitch = 1;
params.hswitch = 1;
params.oswitch = 1;
params.icswitch = 2;

params.index = 0;

inputStruct.params = params;
RunChase1D(inputStruct);
