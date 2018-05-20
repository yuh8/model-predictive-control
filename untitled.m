clearvars
clc
N=20;
M=10;
X=randn(N,M);
Q = X'*X;
f=randn(M,1);
G = randn(5,10);
g = ones(5,1);
for i =1:100000
tic
x = qpFast(Q,f,G,g);
toc;
end

L = chol(Q);
M = G/L;
v = L'\f;
d = g + M*v;

C = -[M';d'];
b = [zeros(10,1);1];

% x = nnls(C,b);


