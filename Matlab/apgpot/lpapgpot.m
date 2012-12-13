function [x,f,R] = lpapgpot(A,b,c,varargin)
%
% lpapgpot sets up the linear program given
% by A,b,c and then solves it using
% the first order solver apgpot
%
%   min_x  c'*x
%   s.t.    A*x = b
%            x >= 0
%
%   4th argument can be a struct "pars" with
%   options and parameters, it may be omitted
%   as all fields have default values set below
%
%    pars.maxit:      Max. number of iterations allowed
%    pars.maxitls:    Max. number of line search backtracks
%    pars.psi:        The function psi(z) ( =g(z) )
%    pars.vart0:      Allow increase in step size (1 or 0)
%    pars.resweight:  Weight on residuals (0, 1 or 2)
%                     0 = no weighting
%                     1 = weight as used in Monteiro/Lan paper
%                     2 = My weighting, which seems to be better
%    pars.chweight:   Allow adaptive change of weights on the fly
%                     0 = no change allowed
%                     N = weights changed every N iterations, N>0
%                         should be on the order 100--1000
%    pars.yfeas:      Keep y feasible (1 or 0)
%    pars.mainstep:   Main step-type (only 'prox' available now)
%    pars.accel:      Use acceleration (1 or 0)
%    pars.echo:       Print-level (1 or 0)
%    pars.prtiter:    Info only printed every pars.prtiter iterations
%    pars.tol:        Tolerance level, relative to intial residuals
%    pars.beta:       Backtracking constant
%    pars.rho:        The parameter in the potential function
%

if nargin < 3
    error('must input at least A, b and c');
end
A     = sparse(A);
b     = sparse(b);
c     = sparse(c);
[m,n] = size(A);
pars  = struct;
K.l   = n;
if nargin > 3 && ~isempty(varargin{1})
    pars = varargin{1};
end
if nargin > 4
    K = varargin{2};
end
if ~isfield(K,'l')
    K.l = n;
end
if ~isfield(K,'f')
    K.f = 0;
end
if n ~= (K.l+K.f)
    error('mismatch size(A,2) vs K.n+K.f');
end


% set defaults:
defaults = {...
    'maxit',100000;...
    'maxitls',100;...
    'psi','psi_log2norm';...
    'inexactprox',0;...
    'vart0',1;...
    'resweight',2;...
    'chweight',0;...
    'yfeas',0;...
    'mainstep','prox';...
    'accel',1;...
    'echo',1;...
    'prtiter',500;...
    'tol',1e-3;...
    'tolnorm',inf;...
    'beta',0.5;...
    'beta3',0.9995;...
    'rho',4*n;... 
    'gamma',1;...
    'hom',0;...
    'stoponlyresx',0;...
    };

for k = 1:size(defaults,1)
    if ~isfield(pars,defaults(k,1))
        pars.(defaults{k,1}) = defaults{k,2};
    end
end

% residual indices
pars.resPidx = 1:m;
pars.resDidx = pars.resPidx(end) + (1:n);
pars.resGidx = pars.resDidx(end) + 1;

% aux
Zmn = sparse(m,n);
Znn = sparse(n,n);
Zmm = sparse(m,m);
In  = speye(n);

% set up M and q from A,b,c:
if ~pars.hom %not homogeneous model
    n2    = K.f;
    n1    = n-n2;
    Af    = A(:,(n1+1):(n1+n2));
    Ap    = A(:,1:n1);
    En1   = speye(n1);
    Zn1n1 = sparse(n1,n1);
    Zn1n2 = sparse(n1,n2);
    Zn2n2 = sparse(n2,n2);
    Zmn1  = sparse(m,n1);    
    
    q     = [b;c;0];
    
    M     = [Ap,Af,Zmn1,Zmm;...
            Zn1n1,Zn1n2,En1,Ap';...
            Zn1n2',Zn2n2,Zn1n2',Af';...
            c',Zn1n1(1,:),-b'];

    pars.bc = [c',Zn1n1(1,:),Zmm(1,:);...
        Znn(1,:),Zn1n1(1,:),b'];

    % n1+n2+n1 first variables must be positive
    pars.pidx = (n1+n2+n1);

    % initial point:
    xy       = ones(n1+n2,1);
    sy       = ones(n1,1);
    yy       = ones(m,1);
    pars.zy0 = [xy;sy;yy];



else % homogeneous model (DOES NOT WORK CORRECTLY RIGHT NOW!)
    
    M = [A,-b,Zmn,Zmn(:,1),Zmm;...
        Znn,-c,In,Znn(:,1),Zmn';...
        -c',0,Znn(1,:),-1,b'];
    q = sparse(m+n+1,1);
    
    pars.bc = [c',0,Znn(1,:),0,Zmm(1,:);...
        Znn(1,:),0,Znn(1,:),0,b'];
    
    pars.pidx     = 2*n + 2;
    pars.tauidx   = n+1;
    pars.kappaidx = 2*n+2;
    
   
    % initial point (just ones)
    tauy   = 1;
    kappay = 1;
    xy       = ones(n,1);
    sy       = xy;
    yy       = ones(m,1);
    pars.zy0 = [xy;tauy;sy;kappay;yy];

end


% norm to use in stopping criteria:
pars.tolnorm = 2;


% residual normalizers
cnorm   = norm(c,pars.tolnorm);
bnorm   = norm(b,pars.tolnorm);
pars.dD = max(cnorm,1);
pars.dP = max(bnorm,1);

% weights on the search direction
% default: no weights:
pars.wP = 1;
pars.wD = 1;
pars.wG = 1;
w = ones(pars.resGidx(end),1);

% if pars.resweight > 0, use weights:
if pars.resweight
    switch pars.resweight
        case 1 % seems to be best for random problems
            pars.wP = 1/(max(1,pars.dP));
            pars.wD = 1/(max(1,pars.dD));
            pars.wG = 1/(max(1,bnorm+cnorm));
        case 2 % seems to be far better for NETLIB problems
            pars.wP = 1/(max(1,norm(A*ones(n,1))));
            pars.wD = 1/(max(1,norm(A'*ones(m,1))));
            pars.wG = 1/(max(1,bnorm+cnorm));
        case 3
            ww = sparse((sum(abs(M),2)).^2);
            pars.wP = ww(pars.resPidx);
            pars.wD = ww(pars.resDidx);
            pars.wG = ww(pars.resGidx);            
    end
    % assign weights to w:
    w(pars.resPidx) = pars.wP;
    w(pars.resDidx) = pars.wD;
    w(pars.resGidx) = pars.wG;
else % if pars.resweight is 0, force no
    % change in the weights (keep them at ones):
    pars.chweight = 0;
end

w = w/pars.gamma;

% run apgpot:
[xx,ff,R] = apgpot(M,q,w,pars);

% output
zsol      = xx;
x         = zsol(1:n);
R.sol.z   = zsol;
R.sol.x   = x;
R.sol.s   = zsol(n+1:(n+n1));
R.sol.y   = zsol((n1+n2+n1):end);
f         = c'*x;
R.M       = M;
R.q       = q;
