function [vd,vq]  = EVAIimplement_Torque_constrained(idset,Tset,id,iq,we,Vdc,dV,Ap,Bp,Cp,Gp,C2,Phi_F,Phi_R,Phi_Gam,H,G,L,M)
xmcur = [id;iq];
Nm = size(Bp,2); % Number of manipulated variable
Ns = size(Ap,1); % Number of state variables
Nout = size(Cp,1); % Number of output variable
Nc = size(C2,2)/Nm; % Control horizon
Nic = Nc/2;
%% Define persistent variables that are local to the function
persistent Xfcur ucur wcur
if isempty(Xfcur)
    Xfcur = zeros(Ns+Nout,1);
end
if isempty(ucur)
    ucur = zeros(Nm,1);
end
if isempty(wcur)
    wcur = 0;
end

%% EVAI implement
yset = [idset;Tset];
QVdq = repmat([0.4;0.9],Nic,1);
QdV = repmat([1.0;1.0],Nic,1);
deltaw = we-wcur;
wcur = we;
%% Solving using QP with constraint
C1 = repmat(eye(Nm),Nic,1);
Umax = QVdq*Vdc;
Umin = -QVdq*Vdc;
dUMax = QdV*dV;
dUMin = -QdV*dV;
fx = Phi_F*Xfcur - Phi_R*yset + Phi_Gam*deltaw;
g = [ -dUMin;
       dUMax;
      C1*ucur - Umin;
     -C1*ucur + Umax];
DeltaU = qpFast(H,fx,G,g,L,M);
%% update
deltau = DeltaU(1:Nm,1);
u = ucur + deltau;
ucur = u;
xm_old = xmcur;
xm = Ap*xmcur + Bp*u + Gp*we;% + Ep;
y = Cp*xm;% + D;
Xfcur = [xm-xm_old;y];
vd = u(1);
vq = u(2);

%% Fast QP
function x = qpFast(Q,c,G,g,L,M)
%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the following quadratic programming problem using
% Non-negative least square active set method
% x= argmin_{x} 0.5 * x'*Q*x + c'*x
%   s.b. Gx <= g
% Copyrigh 2018 Hongyang Yu
%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the size of constraints
q = size(G,1);

% Define alternative formulation as non-negative least square
v = L'\c;
d = g + M*v;

%% Initialize the active set method
P = false(q,1); % Passive set
R = true(q,1); % Active set
qzeros = zeros(q,1);
y = qzeros;
gam = 1;

%% Terminating conditions
epsilon = 1e-3; % termination condition
non_zero = 1e-7; % Approximated non-zero boundary
itmax = 3*q;
in_iter = 0;

coder.varsize('Mp');
coder.varsize('yp');
coder.varsize('dp');
Mp = M;
yp = y;
dp = d;
wy = qzeros;

% Residual for the non-constrained least square problem
w = M*(Mp'*yp) + (gam + dp'*yp)*dp;
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
    A = Mp*Mp' + dp*dp';
    b = -gam*dp;
    S(P) = A\b;
    while any(S(P)<=0) && in_iter<itmax
        in_iter = in_iter+1;
        temp = (S<=0) & P;
        alpha = min(y(temp)./(y(temp)-S(temp)));
        y = y + alpha* (S-y);
        I = P & ((abs(y)<=non_zero));
        gam = gam - norm(d(I),1);
        P(I) = false;
        R = ~P;
        Mp = M(P,:);
        dp = d(P);
        A = Mp*Mp' + dp*dp';
        b = -gam*dp;
        S(P) = A\b;
        S(R) = 0;
    end
    y = S;
    yp = y(P);
    w = M*(Mp'*yp) + (gam + dp'*yp)*d;
end

cond = norm(Mp'*yp)^2 + (gam + dp'*yp)^2;
if cond > non_zero % Approximately nonzero
    lambda = -1/(gam+dp'*yp)*yp;
    u = Mp'*lambda;
    x = L\(u-v);
else
    % Infeasible solution returns 0 vectors
    x = zeros(size(Q,1),1);
end