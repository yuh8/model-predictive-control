clc
clearvars

%% Adjustable parameters
Td = 1e-04; % Discretisation interval of EVAI
Ts = 0.01; % Motor response sampling time
Qid = 0.1;
Qt = 0.5;
Coeff = 0.001;
Np = 5; % Prediction horizon
Nc = 2; % Control horizon

%%
numP = 8; %number of poles
N = 3500; %revolutional speed(r/m) 
Ld = 2.67e-004; %d-axis inductance(H)
Lq = 6.95e-004; %q-axis inductance(H)
Iq0 = 28; %Reference Iq 
Id0 = 0; %Reference Id
fai = 0.09; %flux linkage(Wb)
Vdc = 320; %DC voltage(V)
dV = 100;% Limit voltage increment
R = 0.019; %phase resistance(ohm)
J = 0.032; %inertia(kgm^2)
Dv = 0.001; %viscosity (Nm/(rad/s))
cfreq = 10000; %carrier frequency of PWM(Hz)
wm = N/60*2*pi; %Mechanical Speed(rad/sec)
we0 = wm*numP/2; %Electrrical speed

%%
A = [-R/Ld Lq/Ld*we0;
     -Ld/Lq*we0 -R/Lq];
B = [1/Ld 0;
     0 1/Lq];
C = [1 0;
    3/2*(numP/2)*(Ld-Lq)*Iq0 3/2*(numP/2)*(fai+(Ld-Lq)*Id0)];
Gi = [Lq/Ld*Iq0;
    -Ld/Lq*Id0-fai/Lq];
D = [0;
     -(3/2*(numP/2)*(Ld-Lq)*Id0*Iq0)];
E = [-Lq/Ld*Iq0*we0;
           Ld/Lq*Id0*we0];
Ap= expm(A*Td);
Bp = A\(expm(A*Td)-eye(size(A)))*B;
Cp = C;
Gp = A\(expm(A*Td)-eye(size(A)))*Gi;
Ep = A\(expm(A*Td)-eye(size(A)))*E;
Ap = round(Ap,4);
Bp = round(Bp,4);
Cp = round(Cp,4);
Gp = round(Gp,4);
D  = round(D,4);
Ep = round(Ep,4);

[Phi_Phi, Phi_F, Phi_R, Phi_Gam, C2, H, G, L, M] = MPCmodel_Torque(Ap,Bp,Cp,Gp,Np,Nc,Qid,Qt,Coeff);