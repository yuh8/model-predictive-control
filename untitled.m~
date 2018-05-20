clearvars
clc
N=20;
M=10;
X=randn(N,M);
Q = X'*X;
f=randn(M,1);
G = randn(5,10);
g = ones(5,1);
% x = qpFast(Q,f,G,g);
L = chol(Q);
M = G/L;
v = L'\f;
d = g + M*v;

C = -[M';d'];
b = [zeros(10,1);1];

x = lsqnonneg(C,b);


