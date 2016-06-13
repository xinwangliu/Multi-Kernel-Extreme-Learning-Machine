function [radius, optimalBeta]= minBallRadius(K)

n=size(K,1);
H=2*K;
f=diag(K);
A=ones(n,1);
b=1;
C=1;

[beta, lambda, pos] = monqp(H,f,A,b,C);
optimalBeta=zeros(n,1);
optimalBeta(pos)=beta;

radius = -0.5*optimalBeta'*H*optimalBeta + f'*optimalBeta;
