function [x,f,R] = apgpot(M,q,w,pars)

%
% aims at solving   min_z || Mz - q ||_2^2, s.t. z_I >= 0
% by solving
%  min_z  psi - sum(log(z_I))
%
% where 
%    I is a index set, i.e. just some of the elements of z must be >= 0
%    psi is a some log-potential that is externally specified
%      (see e.g. psi_log2norm)
%    the most obvious example is psi(z) = log( || Mz - q ||_2^2 )
%
% as an example, this can be used to solve LPs 
% by appropriately defining M (see lpapgpot)
%
%
%
% NOTE: In actuality, we are using the potential of the 
%       WEIGHTED 2-norm || W(Mz-q) ||_2^2 as this helps 
%       improve ill-conditioning. See in code exactly how this is 
%       handled. 
%       For LPs, W is just a diagonal matrix with three 
%       diagonal block, one for each residual. 
%       The code now allows W to change adaptively during the iteration
%       depending on the current relative magnitudes of the three
%       residuals. This can improve iteration count substantially
%       but sometimes introduces some instabilities.
%       See the Lan/Monteiro paper for a describtion of using a 
%       static W. 





% for timing:
tic;
stime = toc;

% counters:
NMmul  = 0;
NMtmul = 0;

% set flags:
R.ef     = -9;
R.status = 'maxit';

% init:
zy    = pars.zy0;
zx    = zy;
k     = 0;
t0    = 1;
t1    = 1;
t0m   = 1/0.9;
maxg  = -inf;
PU    = zeros(pars.maxit,3);

% scale q:
% thus in effect letting f(z) = || Mz - gamma*q ||^2
% afterwards, scale z by 1/gamma. This can sometimes 
% be helpful in terms of ill-conditioning
qunscal  = q;
q        = pars.gamma * q;

% printing:
if pars.echo > 0
    fprintf('==================================== APGPOT =======================================\n');
    fprintf('  ITE LSTP F(X)/F0  PROG     PR0/R0   P-RES    D-RES    G-RES    STP      W-ratio \n');
    fprintf('-----------------------------------------------------------------------------------\n');
end

% first residuals:
NMmul  = NMmul + 1;
Mzy    = M*zy;
Mzx    = Mzy;
resy   = Mzy - q;
resx   = resy;
resy0  = Mzy - qunscal;
resx0  = resy0;
resx0n = norm(resx0,2);
dD     = pars.dD;
dP     = pars.dP;
rsP    = max(norm(resx0(pars.resPidx),pars.tolnorm),1);
rsD    = max(norm(resx0(pars.resDidx),pars.tolnorm),1);
rsG    = max(norm(resx0(pars.resGidx),pars.tolnorm),1);

% heuristic for adaptive change of residual weights:
KK    = sqrt(2);
C     = 1 + (KK-1)*(pars.chweight/250)^(3/2);

% first weights:
wP    = pars.wP;
wD    = pars.wD;
wG    = pars.wG;

for j = 1:(pars.maxit+1)

    % NB: remember:
    % only zx stays feasible in cone, i.e. >= 0 for LP
    % zy is allowed to become non-feasible, i.e. <0
    % therefore, only check residuals of zx, NOT zy

    dG  = max( sum(abs(pars.bc*zx))/2 , 1);

    % residuals
    rP  = norm(resx(pars.resPidx),pars.tolnorm);
    rD  = norm(resx(pars.resDidx),pars.tolnorm);
    rG  = norm(resx(pars.resGidx),pars.tolnorm);

    rP  = rP/pars.gamma;
    rD  = rD/pars.gamma;
    rG  = rG/pars.gamma;

    if pars.hom
        rP = rP/zx(pars.tauidx);
        rD = rD/zx(pars.tauidx);
    end

    nrP = rP/rsP;
    nrD = rD/rsD;
    nrG = rG/rsG;

    % measures progress in two different ways:
    % first is relative to initial residuals
    [progr0,i0] = max([nrP,nrD,nrG]);

    % second is like the Lan/Monteiro paper:
    % second should be used when doing comparisons
    % with that paper
    prog        = max([rP/dP,rD/dD,rG/dG]);

    % change residual weights adaptively (if active)
    if pars.chweight > 0
        if mod(j,pars.chweight)==0
            cP = 1/sqrt(C); cD = cP; cG = cP;
            if i0 == 1
                cP = C;
            elseif i0 == 2
                cD = C;
            else
                cG = C;
            end
            wP = cP*wP;
            wD = cD*wD;
            wG = cG*wG;
            w(pars.resPidx,1) = wP;
            w(pars.resDidx,1) = wD;
            w(pars.resGidx,1) = wG;

        end
    end

    % printing:
    resxnorm  = norm(resx,2)/pars.gamma;
    resxnormn = resxnorm/resx0n;
    if pars.echo > 0
        if mod(j-1,pars.prtiter)==0 || pars.echo > 1 || j==1
            maxW = full(max([wP./pars.wP;wD./pars.wD;wG./pars.wG]));
            minW = full(min([wP./pars.wP;wD./pars.wD;wG./pars.wG]));
            fprintf('%5.1i %3.1i %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e  \n',...
                j-1,k,resxnormn,prog,progr0,nrP,nrD,nrG,t1,maxW/minW);
        end
    end

    % record progress, store in PU:
    PU(j,:) = [prog,progr0,resxnormn^2];

    % stopping:
    if pars.stoponlyresx
        if resxnormn < pars.tol
            R.status = 'optimal-resx0';
            R.ef     = 0;
            break;
        end
    else
        if prog < pars.tol
            R.status = 'optimal-prog';
            R.ef     = 0;
            break;
        end
        if progr0 < pars.tol
            R.status = 'optimal-prog0';
            R.ef     = 0;
            break;
        end
    end
    % stopping criteria end

    % safety:
    if isnan(nrP) && isnan(nrD)
        R.status = 'nan';
        R.ef     = -5;
        break;
    end

    % compute prox step:
    resyw        = resy.*w;
    NMtmul       = NMtmul + 1;
    Mtresyw2     = M'*(w.*resyw);
    [psiy,dpsiy] = feval(pars.psi,resyw,Mtresyw2,pars);
    gy           = pars.rho * psiy;
    dgy          = pars.rho * dpsiy;

    % estimate maximal gradient so far:
    maxgrad_candidate = norm(Mtresyw2,'inf');
    if maxgrad_candidate > maxg
        maxg = maxgrad_candidate;
    end

    % store previous x:
    zxprev  = zx;
    Mzxprev = Mzx;

    switch pars.mainstep
        case 'prox'

            % first try step length = t0
            t1 = t0;

            for k = 1:pars.maxitls

                d     = zy - t1*dgy;

                % compute prox(d)
                % (broken into steps to allow viewing
                %  which takes more computational effort):
                tmp   = d.*d;
                tmp   = tmp + 4*t1;
                tmp   = sqrt(tmp);    % sqrt can be slow with extremely 
                tmp   = d + tmp;      % many vareables. In this case, 
                proxd = (1/2)*tmp;    % 3-4 Newton steps to compute the 
                                      % sqrt is faster.
                

                % only apply barrier on pars.pidx first variables:
                % (i.e. those that are required >= 0)
                d(1:pars.pidx) = proxd(1:pars.pidx);

                tmp    = zy - d;
                Gty    = (1/t1)*tmp;
                nGty   = norm(Gty,2)^2;
                gtG    = dgy'*Gty;
                zx     = zy - t1*Gty;

                % if any variable in zx is negative,
                % reduce t1, and skip rest of loop
                if any(zx(1:pars.pidx)<0)
                    t1 = pars.beta * t1;
                    continue;
                end

                NMmul = NMmul + 1;
                Mzx   = M*zx;
                resx  = Mzx - q;
                resxw = resx.*w;
                phix  = feval(pars.psi,resxw,[],pars);
                gx    = pars.rho * phix;

                % When decrease is sufficient, stop. 
                % otherwise, reduce t1
                if gx < gy - t1*gtG + (t1/2)*nGty
                    break;
                else
                    t1 = pars.beta * t1;
                end
            end

            % sometimes it is good to allow trying a
            % slightly larger t0 next time (see TFOCS documentation):
            if pars.vart0
                t0 = t0m*t1;
            else
                t0 = t1;
            end


    end

    % apply acceleration only if pars.accel = 1
    % i.e.: pars.accel=1 turns the algorithm into the
    % accelerated proximal gradient method
    % otherwise it is just the proximal gradient method
    if pars.accel

        % correction search dir:
        Gx   = zx - zxprev;
        MGx  = Mzx - Mzxprev;
        jk   = (j-1)/(j+2);

        % update
        if pars.yfeas
            % if pars.yfeas = 1, then zy is
            % forced to be feasible in each step
            % i.e., zy(pars.pidx) >= 0.
            % this is not strictly necessary, so
            % by default, yfeas = 0
            % again, broken into steps to see computational effort 
            % in each step
            tmp = Gx./zx;
            tmp = -tmp;
            tmp = tmp(1:pars.pidx);
            tmp = max(tmp);
            amx = 1/tmp;
            amx = min(2*abs(amx)/(1+sign(amx)),1);
            t2  = min(jk,pars.beta3*amx);
        else
            t2  = jk;
        end

        % take step
        zy   = zx + t2*Gx;

        % update data (to avoid an extra M or M' multiplication
        Mzy  = Mzx + t2*MGx;
        resy = Mzy - q;
    else
        % if no acceleration, simply update trivially:
        zy   = zx;
        Mzy  = Mzx;
        resy = Mzy - q;
    end

end

% output #iterations:
R.iter = j-1;

PU(j+1:end,:) = [];
R.f           = PU(:,1);
R.f0          = PU(:,2);
R.res         = PU(:,3); %this is f(x)/f(x0)!

% printing:
if pars.echo > 0
    fprintf('%5.1i %3.1i %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e  \n',...
        j-1,k,resxnorm/resx0n,prog,progr0,nrP,nrD,nrG,t1,maxW/minW);
    fprintf('-----------------------------------------------------------------------------------\n');
    fprintf(' %s, #Iters = %i \n #M-mult = %i, #MT-mult = %i, SUM = %i \n',...
        R.status,R.iter,NMmul,NMtmul,NMmul+NMtmul);
    fprintf('===================================================================================\n');
end

R.maxgrad  = maxg;
R.sol.agpn = zx/pars.gamma;
x          = R.sol.agpn;
f          = resxnorm^2;
R.ttime    = toc - stime;

