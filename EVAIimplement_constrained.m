function [vd,vq]  = EVAIimplement_constrained(idset,wset,id,iq,we,Vdc,C2,Ap,Bp,Cp,Phi_Phi,Phi_F,Phi_R,Coeff)
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
Umax = C1*Vdc/sqrt(3);
Umin = -C1*Vdc/sqrt(3);
H = Phi_Phi + Coeff*eye(size(Phi_Phi));
fx = Phi_F*Xfcur - Phi_R*yset;
bmin = Umin - C1*ucur;
bmax = Umax - C1*ucur;
[DeltaU,~,~,~,~] = qp_constrained(H,fx,C2,bmin,bmax);
%% update
deltau = DeltaU(1:Nm,1);
u = ucur + deltau;
ucur = u;
xm_old = xmcur;
xm = Ap*xmcur + Bp*u;
y = Cp*xm;
Xfcur = [xm-xm_old;y];
vd = u(1);
vq = u(2);
