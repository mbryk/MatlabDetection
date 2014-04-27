function [rate,PF,PD] = DetectorF(threshold,A,sigma,sigma2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

iter = 500;
Ri=zeros(iter,1);PFi=Ri;PDi=Ri;

for i=1:iter
% Create Target
target = randi(10,100,1);
target(target<=8)=0;
target(target>8)=1;

% Create Signal
y = target;
y(target==0) = normrnd(A,sigma2,size(y(y==0))); % Y=A+Z
y(target==1) = normrnd(A,sigma,size(y(y==1))); % Y=A+X

% MAP Detection 
y = -((y-A).^2);
d=y;
d(y<threshold) = 0;
d(y>=threshold) = 1;

Ri(i) = sum(target~=d)/length(target);
H0 = d(target==0);
H1 = d(target==1);
PFi(i) = sum(H0)/length(H0);
PDi(i) = sum(H1)/length(H1);
end

rate = mean(Ri);
PF = mean(PFi);
PD = mean(PDi);
end