function model = lds_EM(U,Y,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%   U: input signal
%   Y: output signal
%   M: dimension of Latent space
%   iter: Number of EM iterations
%Output:
%   Estimated LDS model parameters
% Written and tested by Hongyang Yu @ yhytxwd@gmail.com

if nargin<2
    error('Not enough number of input argument')
end

if nargin<3
    M = 2;
    iter = 100;
elseif nargin<4
    M = varargin{1};
    iter = 100;
else
    M = varargin{1};
    iter = varargin{2};
end

%% Determine the dimensions of input,output and state
[O,T] = size(U);
[D,~] = size(Y);

%% Parameter initialization
A = 0.1*eye(M);
B = 0.1*ones(M,O);
C = 0.1*ones(D,M);
Q = 0.001*eye(M,M);
R = 0.001*eye(D,D);
mu = zeros(M,1);
v = 0.01*eye(M);

LL = zeros(1,iter);
%% EM iteration
for ii = 1:iter
    % Forward recursion
    [Xcur, Vcur, Vpre,LL(ii)] = KalmanFilter(U,Y,A,B,C,mu,v,Q,R);

    % Backward recursion and E-step
    [A_sum,B_sum,C_sum,R_sum,U_sum,Pt_sum,Pt_sum1,Pt_sum2,Ptt1_sum] = KalmanSmoother(U,Y,Xcur,Vcur,Vpre,A,B,Q);
    
    % M-step
    C = C_sum/Pt_sum;
    R = diag(diag((R_sum - C*C_sum')))/T;
    A = (Ptt1_sum - B*A_sum)/Pt_sum1;
    B = (B_sum - A*A_sum')/U_sum;
    Q = (Pt_sum2 - Ptt1_sum*A' - A*Ptt1_sum' + A*Pt_sum1*A'...
        - B_sum*B' - B*B_sum' +A*A_sum'*B' + B*A_sum*A' + B*U_sum*B')/(T-1);
    mu = Xcur(:,1);
    v = Vcur(:,:,1);
end
model.A = A;
model.B = B;
model.C = C;
model.Q = Q;
model.R = R;
model.mu = mu;
model.v = v;
model.LL = LL;


