function [A_sum,B_sum,C_sum,R_sum,U_sum,Pt_sum,Pt_sum1,Pt_sum2,Ptt1_sum] = KalmanSmoother(U,Y,Xcur,Vcur,Vpre,A,B,Q)
%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%   Xcur is the filtered state vector
%   Vcur is the variance of the estimated state vector
%   Vpre is the predicted variance of the state vector
%Output:
%   Summation terms for M-step
% Writen and tested by Hongyang Yu @yhytxwd@gmail.com

% Determine the size of matrix
[M,T] = size(Xcur);
[O,~] = size(U);

% Memory pre-allocation for speed
Xfin = zeros(M,T);
Vfin = zeros(M,M,T);
Ptt1_sum = zeros(M,M);
Pt_sum1 = zeros(M,M);
Pt_sum2 = zeros(M,M);
A_sum = zeros(O,M);
B_sum = zeros(M,O);
Q_sum = zeros(M,M);
U_sum = zeros(O,O);

% Initialize the smoother
Xfin(:,T) = Xcur(:,T);
Vfin(:,:,T) = Vcur(:,:,T);

Pt_sum = Vcur(:,:,T) + Xcur(:,T)*Xcur(:,T)';
C_sum = Y(:,T)*Xfin(:,T)';
R_sum = Y(:,T)*Y(:,T)';

% Kalman smoother
for t = T:-1:2
    % Expectation term for E-step
    J = Vcur(:,:,t-1)*A'/Vpre(:,:,t);
    S = Vcur(:,:,t-1) - J*Vpre(:,:,t)*J';
    Xfin(:,t-1) = Xcur(:,t-1) + J*(Xfin(:,t) - A*Xcur(:,t-1)) - S*A'*inv(Q)*B*U(:,t-1);
% Xfin(:,t-1) = Xcur(:,t-1) + J*(Xfin(:,t) - A*Xcur(:,t-1));
    Vfin(:,:,t-1) = Vcur(:,:,t-1) + J*(Vfin(:,:,t)-Vpre(:,:,t))*J';
    Vtt1 = Vfin(:,:,t)*J';
    
    % Summation terms for M-step
    Pt_sum = Pt_sum + Vfin(:,:,t-1) + Xfin(:,t-1)*Xfin(:,t-1)';
    Pt_sum1 = Pt_sum1 + Vfin(:,:,t-1) + Xfin(:,t-1)*Xfin(:,t-1)';
    Pt_sum2 = Pt_sum2 + Vfin(:,:,t) + Xfin(:,t)*Xfin(:,t)';
    Ptt1_sum = Ptt1_sum + Vtt1 + Xfin(:,t)*Xfin(:,t-1)';
    A_sum = A_sum + U(:,t-1)*Xfin(:,t-1)';
    B_sum = B_sum + Xfin(:,t)*U(:,t-1)';
    C_sum = C_sum + Y(:,t-1)* Xfin(:,t-1)';
    R_sum = R_sum + Y(:,t-1)*Y(:,t-1)';
    U_sum = U_sum + U(:,t-1)*U(:,t-1)';
end