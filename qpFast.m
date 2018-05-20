function x = qpFast(Q,c,G,g)
%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the following quadratic programming problem using
% Non-negative least square active set method
% x= argmin_{x} 0.5 * x'*Q*x + c'*x
%   s.b. Gx <= g
% Copyrigh 2018 Hongyang Yu
%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the size of constraints
q = size(G,1);

%% Offline computation which can be taken out of this function!!!
% Define alternative formulation as non-negative least square
L = chol(Q);
M = G/L;
v = L'\c;
d = g + M*v;

%% Initialize the active set method
P = false(q,1); % Passive set
R = true(q,1); % Active set
qzeros = zeros(q,1);
y = qzeros;
gam = 1;

%% Terminating conditions
epsilon = 1e-4; % termination condition
non_zero = 1e-6; % Approximated non-zero boundary
itmax = 100;

Mp = M;
yp = y;
dp = d;
wy = qzeros;

% Residual for the non-constrained least square problem
w = M*(M'*y) + (gam + d'*y)*d;
while any(w<-(gam + dp'*yp)*epsilon) ...
        && any(R)...
        && (norm(Mp'*yp)^2 + (gam + dp'*yp)^2>non_zero)
    
    S = qzeros;
    wy(P) = inf;
    wy(R) = w(R);
    [~,i] = min(wy);
    P(i) = true;
    R(i) = false;
    gam = gam + abs(d(i));
    
    Mp = M(P,:);
    dp = d(P);
    yp = y(P);
    A = Mp*Mp' + dp*dp';
    b = -gam*dp;
    S(P) = A\b;
    in_iter = 0;
    while any(S(P)<=0) || in_iter<itmax
        in_iter = in_iter+1;
        temp = (S<=0) & P;
        alpha = min(y(temp)./(y(temp)-S(temp)));
        y = y + alpha* (S-y);
        I = P & ((abs(y)<=non_zero));
        gam = gam - norm(d(I),1);
        R = I | R;
        P = ~R;
        Mp = M(P,:);
        dp = d(P);
        yp = y(P);
        % The main reason this alg is efficient as only P terms are stored
        % everytime
        A = Mp*Mp' + dp*dp';
        b = -gam*dp;
        S(P) = A\b;
        S(R) = 0;
    end
    y = S;
    w = M*(Mp'*yp) + (gam + dp'*yp)*d;
end

cond = norm(Mp'*yp)^2 + norm(gam + dp'*yp)^2;
if cond > non_zero % Approximately nonzero
    lambda = -1/(gam+dp'*yp)*yp;
    u = Mp'*lambda;
    x = L\(u-v);
else
    % Infeasible solution returns 0 vectors
    x = zeros(size(Q,1),1);
end



