function [psi,dpsi] = psi_log2norm(res,Mtres,pars)

f    = norm(res,2)^2;
psi  = log(f);
dpsi = (2/f)*Mtres;
