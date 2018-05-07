function [u, y, xm,Xf] = MPCimplement(ucur,yset,Xfcur,xmcur,model,Coeff,Phi_Phi,Phi_F,Phi_R)
%%Calculating cost function
Ni = length(ucur);
Ad = model.Ad;
Bd = model.Bd;
Cd = model.Cd;

temp = Phi_Phi+Coeff*eye(size(Phi_Phi));
invtemp = choleskyINV(temp);
DeltaU = invtemp*(Phi_R*yset-Phi_F*Xfcur);
deltau = DeltaU(1:Ni,1);
u = ucur + deltau;
xm_old = xmcur;
xm = Ad*xmcur+Bd*u;
y = Cd*xm;
Xf = [xm-xm_old;y];

% Stable inversion of massive positive definite matrix
function invX = choleskyINV(X)
L = chol(X,'lower');
invX = L'\(L\eye(size(X)));

