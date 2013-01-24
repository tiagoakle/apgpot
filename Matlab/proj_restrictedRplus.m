function op = proj_restrictedRplus( n,m )

  %TFOCS compatible projection
  %Returns an implementation of the indicator function for 
  %the non negative orhtant of the entries from 1:n and from n+m+1:2n+m
  %for a variable x of size 2*n+m

if nargin ~= 2,
  error('Requires m,n as arguments');
elseif ~isa( m, 'double' ) || numel( m ) ~= 1 || ~isreal( m ),
    error( 'The argument m must be a real scalar.' );
elseif ~isa( n, 'double' ) || numel( n ) ~= 1 || ~isreal( n ),
    error( 'The argument n must be a real scalar.' );
else
  op = @proj_Rplus_impl;
end

function [ v, x ] = proj_Rplus_impl( x, t )
v = 0;
switch nargin,
	case 1,
		if nargout == 2,
			error( 'This function is not differentiable.' );
        elseif any( x(1:n) < 0 )||any(x(n+m+1:2*n+m) < 0),
            v = Inf;
        end
	case 2,
		x(1:n)         = max( x(1:n), 0 );
        x(n+m+1:2*n+m) = max( x(n+m+1:2*n+m), 0 );
	otherwise,
		error( 'Not enough arguments.' );
end
end
end
% TFOCS v1.2 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
