%% MATLAB Detection Exercise
%  Mark Bryk and Yaron Tokayer
%  ECE 302 - Stochastics and Probability
%  5/1/14
%

%% 
clc, clear, close all

%% Part 2 - Photon Detection


P0 = .8; P1 = 1-P0;
P = P0*10;
bits = randi(10,1000,1);
bits(bits<=P)=0;
bits(bits>P)=1;

lambda0=50; lambda1=500;
mu0 = 1/lambda0; mu1 = 1/lambda1;
y = bits;
y(bits==0) = exprnd(mu0,size(y(bits==0)));
y(bits==1) = exprnd(mu1,size(y(bits==1)));

C = [0 1; 1 0];
eta = (C(2,1)-C(1,1))/(C(1,2)-C(2,2)) * (P0/P1);
threshold = log(lambda1/(eta*lambda0))/(lambda0+lambda1);
d=y;
d(y>threshold) = 0;
d(y<=threshold) = 1;

rate = sum(d~=bits)/length(bits);