function [Phi_Phi, Phi_F, Phi_R, Phi_Gam, C2, H, G, L, M] = MPCmodel_Torque(Ad,Bd,Cd,Gd,Np,Nc,idw,Tw,Coeff)

[Nout,~] = size(Cd); % Number of output variables
[Ns,~] = size(Ad); % Number of state variables
[~,Nm] = size(Bd); % Number of manipulated variables

%%SSM Augmentation
A = [Ad, zeros(Ns,Nout);Cd*Ad, eye(Nout,Nout)];
B = [Bd;Cd*Bd];
G = [Gd;Cd*Gd];
C = [zeros(Nout,Ns), eye(Nout,Nout)];

%%Calculating key matrics in cost function J
F = C*A;
Temp_F = F;
Phi = C*B;
Gam = C*G;
temp = C;
Q = diag(repmat([idw,Tw],1,Np));
for ii=1:Np-1
    temp = temp*A;
    Phi = [Phi; temp*B];
    Gam = [Gam; temp*G];
    Temp_F=Temp_F*A;
    F= [F;Temp_F];
end

temp = Phi;
[s1,s2] = size(C*B);
temp1 = repmat(eye(Nm),Nc,1);
C2 = temp1;
for ii = 1:Nc-1
    temp2 = [repmat(zeros(s1,s2),ii,1); temp(1:end-s1*ii,:)];
    Phi = [Phi,temp2];
    temp3 = [repmat(zeros(Nm,Nm),ii,1); temp1(1:end-Nm*ii,:)];
    C2 = [C2,temp3];
end

Phi_Phi= Phi'*Q*Phi;
Phi_F= Phi'*Q*F;
BarRs=repmat(eye(Nout),Np,1);
Phi_R=Phi'*Q*BarRs;
Phi_Gam = Phi'*Q*Gam;

H = Phi_Phi + Coeff*eye(size(Phi_Phi));
C2 = C2(1:(Nc*Nm/2),:);
G = [-C2;C2;-C2;C2];
L = chol(H);
M = G/L;
