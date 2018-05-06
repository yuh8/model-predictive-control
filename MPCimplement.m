function [u, y, xm,Xf] = MPCimplement(ucur,yset,Xfcur,xmcur,ssmodel,Coeff,Phi_Phi,Phi_F,Phi_R)
%%Calculating cost function
Ni = length(ucur);
Ap = ssmodel.A;
Bp = ssmodel.B;
Cp = ssmodel.C;

DeltaU = inv(Phi_Phi+Coeff*eye(size(Phi_Phi)))*(Phi_R*yset-Phi_F*Xfcur);
deltau = DeltaU(1:Ni,1);
u = ucur + deltau;
xm_old = xmcur;
xm = Ap*xmcur+Bp*u;
y = Cp*xm;
Xf = [xm-xm_old;y];
