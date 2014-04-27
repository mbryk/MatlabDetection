function [ rate,PF,PD ] = PhotonDetector(threshold,lambda0,lambda1 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

iter = 500;
Ri=zeros(iter,1);PFi=Ri;PDi=Ri;

for i=1:iter

P0 = .8; P1 = 1-P0;
P = P0*10;
bits = randi(10,1000,1);
bits(bits<=P)=0;
bits(bits>P)=1;

mu0 = 1/lambda0; mu1 = 1/lambda1;
y = bits;
y(bits==0) = exprnd(mu0,size(y(bits==0)));
y(bits==1) = exprnd(mu1,size(y(bits==1)));

d=y;
d(y>threshold) = 0;
d(y<=threshold) = 1;

rate = sum(d~=bits)/length(bits);

end

