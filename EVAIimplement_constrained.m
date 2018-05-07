function [vd,vq]  = EVAIimplement_constrained(idset,wset,id,iq,we,V_DC,C2,Ap,Bp,Cp,Phi_Phi,Phi_F,Phi_R,Coeff)
xmcur = [id;iq;we];
Nm = size(Bp,2); % Number of manipulated variable
Ns = size(Ap,1); % Number of state variables
Nout = size(Cp,1); % Number of output variable
Nc = size(C2,1)/Nm; % Control horizon
%% Define persistent variables that are local to the function
persistent Xfcur ucur
if isempty(Xfcur)
    Xfcur = zeros(Ns+Nout,1);
end
if isempty(ucur)
    ucur = zeros(Nm,1);
end

%% MPC implement
yset = [idset;wset];
%% Solving using QP with constraint
C1 = repmat(eye(Nm),Nc,1);
Umax = C1*V_DC/sqrt(3);
Umin = -C1*V_DC/sqrt(3);
H = Phi_Phi+Coeff*eye(size(Phi_Phi));
fx = -2*(Phi_R*yset-Phi_F*Xfcur);
A = C2;
bmin = Umin-C1*ucur;
bmax = Umax-C1*ucur;
[DeltaU,~,~,~,~] = qp_constrained(H,fx,A,bmin,bmax);
%% update
deltau = DeltaU(1:Nm,1);
u = ucur + deltau;
ucur = u;
xm_old = xmcur;
xm = Ap*xmcur+Bp*u;
y = Cp*xm;
Xfcur = [xm-xm_old;y];
vd = u(1);
vq = u(2);
