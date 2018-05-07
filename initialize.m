clc
clearvars

%% Adjustable parameters
Td = 1e-4;
Ts = 1e-4; % sampling time of PMSM model and MPC discritization
Qid = 50;
Qw = 0.5;
Coeff = 0.05;
Np = 20; % Prediction horizon
Nc = 10; % Control horizon

%%
numP = 8; %number of poles
N = 3000; %revolutional speed(r/m) 
Ld = 2.67e-004; %d-axis inductance(H)
Lq = 6.95e-004; %q-axis inductance(H)
Iq0 = 12.275; %Reference Iq 
Id0 = -1.725; %Reference Id
fai = 0.09; %flux linkage(Wb)
Vdc = 320; %DC voltage(V)
R = 0.019; %phase resistance(ohm)
J = 0.02; %inertia(kgm^2)
Vs = 0.006; %viscosity (Nm/(rad/s))
cfreq = 10000; %carrier frequency of PWM(Hz)
Temp = 293.15; %Temperature(K)
Tm = 0; %Load(Nm)
wm = N/60*2*pi; %frequency of speed control(rad/sec)
we0 = wm*numP/2; %frequency of current control(rad/sec)
Kt = 0.629; %torque constant of motor(Nm/A)

%%
A = [-R/Ld Lq/Ld*we0 Lq/Ld*Iq0;
     -Ld/Lq*we0 -R/Lq -(Ld/Lq*Id0+fai/Lq);
     0 3*(numP/2)^2*fai/(2*J) -Vs/J];
B = [1/Ld 0;
     0 1/Lq;
     0 0];
C = [1 0 0;
     0 0 1];
D = [0 0;0 0];

% Used analytical solution rather than matlab discretization 
Ap = expm(A*Ts);
Bp = A\(expm(A*Ts)-eye(size(A)))*B;
Cp = C;
[Phi_Phi, Phi_F, Phi_R] = MPCmodel(Ap,Bp,Cp,Np,Nc,Qid,Qw);
