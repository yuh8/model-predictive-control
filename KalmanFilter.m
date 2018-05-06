function [Xcur, Vcur, Vpre,LL] = KalmanFilter(U,Y,A,B,C,mu,v,Q,R)
%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%   U is the DEMEANED input signal matrix
%   Y is the DEMEANED output data matrix
%   A is the state tranistion matrix
%   B is the input matrix
%   C is the emission matrix
%   mu is the mean of the initial state vector
%   v is the variance of the initial state vector
%   Q is the covariance of state transition noise
%   R is the measurement nosie
%Output:
%   Xcur is the estimated mean of the state vector
%   Vcur is the estimated covariance of the state vector
%   Vpre is the one-step-ahead prediection of variance
%   LL is the sum of log-likelihood
% Written and tested by Hongyang Yu @yhytxwd@gmail.com

% Determine the dimensions of variables
M = length(mu);
[~,T] = size(U);

% Memory pre-allocation for acceleration
Vcur = zeros(M,M,T);
Vpre = zeros(M,M,T);
Xcur = zeros(M,T);
LL = 0;

% Initializing the filter
Vpre(:,:,1) = v;
Xpre = mu;

% Kalman filter
for t = 1:T
    K = Vpre(:,:,t)*C'/(R+C*Vpre(:,:,t)*C');
    Vcur(:,:,t) = Vpre(:,:,t) - K*C*Vpre(:,:,t);
    Xcur(:,t) = Xpre + K*(Y(:,t)-C*Xpre);
    Sigma = C*Vpre(:,:,t)*C' + R;
    LL = LL + loggausspdf(Y(:,t), C*Xpre, Sigma);
    Xpre = A*Xcur(:,t) + B*U(:,t);
    if t<T
        Vpre(:,:,t+1) = A*Vcur(:,:,t)*A' + Q;
    end
end