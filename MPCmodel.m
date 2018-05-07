function [Phi_Phi, Phi_F, Phi_R] = MPCmodel(Ad,Bd,Cd,Np,Nc,idw,ww)

No = size(Cd,1); % Number of output variables
Ns = size(Ad,1); % Number of state variables

%%SSM Augmentation
A = [Ad, zeros(Ns,No);Cd*Ad, eye(No,No)];
B = [Bd;Cd*Bd];
C = [zeros(No,Ns), eye(No,No)];

%%Calculating key matrics in cost function J
F = C*A;
Temp_F = F;
Phi = C*B;
temp = C;
Q = diag(repmat([idw,ww],1,Np));
for ii=1:Np-1
    temp = temp*A;
    Phi = [Phi; temp*B];
    Temp_F=Temp_F*A;
    F= [F;Temp_F];
end

temp = Phi;
[s1,s2] = size(C*B);
for ii = 1:Nc-1
    temp1 = [repmat(zeros(s1,s2),ii,1); temp(1:end-s1*ii,:)];
    Phi = [Phi,temp1];
end

Phi_Phi= Phi'*Q*Phi;
Phi_F= Phi'*Q*F;
BarRs=repmat(eye(No),Np,1);
Phi_R=Phi'*Q*BarRs;

Phi_Phi = round(Phi_Phi,4);
Phi_F = round(Phi_F,4);
Phi_R = round(Phi_R,4);
