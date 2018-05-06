function [Ypre,Rpre] = KalmanPredictor(U,model)
%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%   U is the input signal matrix
%   SS model parameters
%Output:
%   Ypre is the predicted mean of the output
%   Rpre is the predicted variance

% Written and tested by Hongyang Yu @yhytxwd@gmail.com

% Load model paramters
A = model.A;
B = model.B;
C = model.C;
Q = model.Q;
R = model.R;
mu = model.mu;
v = model.v;
% Determine the dimensions of variables
M = length(mu);
[~,T] = size(U);
[D,~] = size(C);

% Memory pre-allocation for acceleration
Vcur = zeros(M,M,T);
Vpre = zeros(M,M,T);
Ypre = zeros(D,T);
Rpre = zeros(D,D,T);

% Initializing the filter
Vpre(:,:,1) = v;
Xpre = mu;

% Kalman predictor
for t = 1:T
    Rpre(:,:,t) = R+C*Vpre(:,:,t)*C';
    K = Vpre(:,:,t)*C'/Rpre(:,:,t);
    Vcur(:,:,t) = Vpre(:,:,t) - K*C*Vpre(:,:,t);
    Ypre(:,t) = C*Xpre;
    Xpre = A*Xpre + B*U(:,t);
    if t<T
        Vpre(:,:,t+1) = A*Vcur(:,:,t)*A' + Q;
    end
end