P = UFget(637);
A = P.A;
b = P.b;
c = P.aux.c;

x0     = [];
z0     = [];
mu     = 1;
opts   = struct;
opts.saveHist = true;
opts.alg      = 'AT'
opts.maxIter  = 1000;
[sol,res] = solver_sLP( c, A,-b, mu);
plot(semilogy(-res.f));
