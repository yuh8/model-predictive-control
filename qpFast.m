function x = qpFast(Q,c,G,g)
%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the following quadratic programming problem using
% Non-negative least square active set method
% J=0.5 * x'*Q*x + c'*x
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
P = []; % Passive set
R = 1:q; % Active set
set = 1:q; % Full set

y = zeros(q,1);
gam = 1;
epsilon = 1e-4; % termination condition
non_zero = 1e-6; % Approximated non-zero boundary

Mp = M;
yp = y;
dp = d;
S = zeros(q,1);
% Residual for the non-constrained least square problem
w = M*(M'*y) + (gam + d'*y)*d;
while any(w<-(gam + dp'*yp)*epsilon) ...
        || length(intersect(P,set))<q || ~isempty(R)...
        || (norm(Mp'*yp)^2 + norm(gam + dp'*yp)^2>non_zero)
    
    [~,i] = min(w(R));
    % Making sure that P and R are always a set rather than a vector
    temp = i;
    i = R(i); R(temp) = [];
    P = [P,i]; P = sort(unique(P));
    gam = gam + abs(d(i));
    
    Mp = M(P,:);
    dp = d(P);
    yp = y(P);
    A = Mp*Mp' + dp*dp';
    b = -gam*dp;
    S(P) = A\b;
    S(R) = 0;
    if any(S(P)<0)
        temp = P(S(P)<=0);
        [~,idx] = min(y(temp)./(y(temp)-S(y(temp))));
        j = temp(idx);
        y = y + y(j) / (y(j) - S(j)) * (s-y);
        [I,idxP,~] = intersect(P,find(y<non_zero));
        P(idxP) = [];
        gam = gam - sum(abs(d(I)));
        Mp = M(P,:);
        dp = d(P);
        yp = y(P);
        % The main reason this alg is efficient as only P terms are stored
        % everytime
        A = Mp*Mp' + dp*dp';
        b = -gam*dp;
        S(P) = A\b;
        R = setdiff(set,P);
        S(R) = 0;
    else
        y = S;
        w= M*(Mp'*yp) + (gam + dp'*yp)*d;
    end
end

cond = norm(Mp'*S(P)')^2 + norm(gam + dp'*S(P))^2;
if cond > non_zero % Approximately nonzero
    lambda = -1/(gam+dp'*yp)*y;
    u = Mp'*lambda(P);
    x = L\(u-v);
else
    % Infeasible solution returns 0 vectors
    x = zeros(size(Q,1),1);
end



