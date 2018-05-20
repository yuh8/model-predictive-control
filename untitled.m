clearvars
clc

tic;
for i =1:100000
N=20;
M=10;
X=randn(N,M);
Q = X'*X;
f=randn(M,1);
G = randn(5,10);
g = ones(5,1);
x = qpFast(Q,f,G,g);
end
toc;
L = chol(Q);
M = G/L;
v = L'\f;
d = g + M*v;

C = -[M';d'];
b = [zeros(10,1);1];

% x = nnls(C,b);


