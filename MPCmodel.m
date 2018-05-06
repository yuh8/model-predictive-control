function [Phi_Phi, Phi_F, Phi_R] = MPCmodel(ssmodel,Np,Nc)
Ap = ssmodel.A;
Bp = ssmodel.B;
Cp = ssmodel.C;

[D,~] = size(Cp);
Ni = D;
[M,~] = size(Bp);

%%SSM Augmentation
A = [Ap, zeros(M,D);Cp*Ap, eye(D,D)];
B = [Bp;Cp*Bp];
C = [zeros(D,M), eye(D,D)];

%%Calculating key matrics in cost function J
F = C*A;
Temp_F = F;
Phi = C*B;
temp = C;

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

Phi_Phi= Phi'*Phi;
Phi_F= Phi'*F;
BarRs=repmat(eye(Ni),Np,1);
Phi_R=Phi'*BarRs;
