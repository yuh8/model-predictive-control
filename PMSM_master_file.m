clc
close all
clearvars

%% Training SSM
% Load training data
% U(:,1) = xlsread('Voltage.xlsx');
U(:,1) = xlsread('Current.xlsx');
Y = xlsread('Speed.xlsx');
U = U';
Y = Y';
Utrain = U(:,1:200);
Ytrain = Y(:,1:200);
Utest = U;
Ytest = Y;

ssmodel = lds_EM(Utrain,Ytrain,1,200);
M = size(ssmodel.A,1);
D = size(ssmodel.C,1);
O = size(ssmodel.B,2);
[Ypre,Rpre] = KalmanPredictor(U,ssmodel);

Tt = 800;
figure
plot(1:Tt,Ytest(1,1:Tt),'b');
hold on
plot(1:Tt,Ypre(1,1:Tt),'r');
ylabel('Speed')
legend('Actual test data','Predicted by our model')

%% MPC implementation
Np = 20;
Nc = 8;
[Phi_Phi, Phi_F, Phi_R] = MPCmodel(ssmodel,Np,Nc);

xm = zeros(M,1);
Xf = zeros(M+D,1);
yset = xlsread('Speed.xlsx',1);

% Set initial control input
u = 10; % u(k-1) =0
u1 = u';
y = 500;
y1 = y';

Coeff = 0.1;
for kk = 1:Tt
    [u, y, xm,Xf] = MPCimplement(u,yset(kk,1)',Xf,xm,ssmodel,Coeff,Phi_Phi,Phi_F,Phi_R);
    f = @() MPCimplement(u,yset(kk,1)',Xf,xm,ssmodel,Coeff,Phi_Phi,Phi_F,Phi_R);
     t =timeit(f)
    u1 = [u1;u'];
    y1 = [y1;y'];
end

figure
subplot(211)
plot(1:Tt, U(1,1:Tt),'b');
hold on
plot(1:Tt, u1(1:Tt,1),'r');
title('Input comparision')
legend('PID', 'MPC')

% subplot(312)
% plot(1:Tt, U(2,1:Tt),'b');
% hold on
% plot(1:Tt, u1(1:Tt,2),'r');
% title('Input comparision')
% legend('PID', 'MPC')

subplot(212)
plot(1:Tt, Y(1:Tt),'b');
hold on
plot(1:Tt, y1(1:Tt),'r');
title('Output comparision')
legend('PID', 'MPC')

