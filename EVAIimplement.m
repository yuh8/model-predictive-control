function [vd,vq]  = EVAIimplement(idset,wset,id,iq,we,Ap,Bp,Cp,Phi_Phi,Phi_F,Phi_R,Coeff)
xmcur = [id;iq;we];
Nm = size(Bp,2); % Number of manipulated variable
Ns = size(Ap,1); % Number of state variables
No = size(Cp,1); % Number of output variable

%% Define persistent variables that are local to the function
persistent Xfcur ucur
if isempty(Xfcur)
    Xfcur = zeros(Ns+No,1);
end
if isempty(ucur)
    ucur = zeros(Nm,1);
end

%% MPC implement
yset = [idset;wset];
% Stable inversion of massive positive definite matrix
temp = Phi_Phi+Coeff*eye(size(Phi_Phi));
L = chol(temp,'lower');
invtemp = L'\(L\eye(size(temp)));
% Compute DeltaU
DeltaU = invtemp*(Phi_R*yset-Phi_F*Xfcur);
deltau = DeltaU(1:Nm,1);
u = ucur + deltau;
ucur = u;
xm_old = xmcur;
xm = Ap*xmcur+Bp*u;
y = Cp*xm;
Xfcur = [xm-xm_old;y];
vd = u(1);
vq = u(2);
