function [x,resnorm,w,outeriter] = nnls(C,d)
n = size(C,2);
% Initialize vector of n zeros and Infs (to be used later)
nZeros = zeros(n,1);
wz = nZeros;

% Initialize set of non-active columns to null
P = false(n,1);
% Initialize set of active columns to all and the initial point to zeros
Z = true(n,1);
x = nZeros;

resid = d - C*x;
w = C'*resid;

% Set up iteration criterion
outeriter = 0;
itmax = 3*n;
tol = 1e-7;
% Outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(Z) > tol)
   outeriter = outeriter + 1;
   % Reset intermediate solution z
   z = nZeros; 
   % Create wz, a Lagrange multiplier vector of variables in the zero set.
   % wz must have the same size as w to preserve the correct indices, so
   % set multipliers to -Inf for variables outside of the zero set.
   wz(P) = -Inf;
   wz(Z) = w(Z);
   % Find variable with largest Lagrange multiplier
   [~,t] = max(wz);
   % Move variable t from zero set to positive set
   P(t) = true;
   Z(t) = false;
   % Compute intermediate solution using only variables in positive set
   z(P) = C(:,P)\d;
   % inner loop to remove elements from the positive set which no longer belong
   iter = 0;
   while any(z(P) <= 0) || iter<itmax
       iter = iter + 1;
       % Find indices where intermediate solution z is approximately negative
       Q = (z <= 0) & P;
       % Choose new x subject to keeping new x nonnegative
       alpha = min(x(Q)./(x(Q) - z(Q)));
       x = x + alpha*(z - x);
       % Reset Z and P given intermediate values of x
       Z = ((abs(x) < tol) & P) | Z;
       P = ~Z;
       z = nZeros;           % Reset z
       z(P) = C(:,P)\d;      % Re-solve for z
   end
   x = z;
   resid = d - C*x;
   w = C'*resid;
end
resnorm = resid'*resid;
