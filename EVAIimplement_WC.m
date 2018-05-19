function [vd,vq]  = EVAIimplement_WC(idset,tset,id,iq,we,Vdc,C2,Ap,Bp,Cp,Gp,D, Ep,Phi_Phi,Phi_F,Phi_R,Phi_Gam,Coeff)
xmcur = [id;iq];
Nm = size(Bp,2); % Number of manipulated variable
Ns = size(Ap,1); % Number of state variables
Nout = size(Cp,1); % Number of output variable
Nc = size(C2,1)/Nm; % Control horizon

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

%% MPC implement
yset = [idset;tset];
% Solving using QP with constraint
deltaw = we - wcur;
wcur = we;
C1 = repmat(eye(Nm),Nc,1);
Umax = ones(Nc*2,1)*Vdc/sqrt(3);
Umin = -ones(Nc*2,1)*Vdc/sqrt(3);
H = Phi_Phi + Coeff*eye(size(Phi_Phi));
fx = Phi_F*Xfcur - Phi_R*yset + Phi_Gam*deltaw;
bmin = Umin - C1*ucur;
bmax = Umax - C1*ucur;
[DeltaU,~,~,~,~] = qp_constrained(H,fx,C2,bmin,bmax);

%% update
deltau = DeltaU(1:Nm,1);
u = ucur + deltau;
ucur = u;
xm_old = xmcur;
xm = Ap*xmcur + Bp*u + Gp*we + Ep;
y = Cp*xm + D;
Xfcur = [xm-xm_old;y];
vd = u(1);
vq = u(2);
