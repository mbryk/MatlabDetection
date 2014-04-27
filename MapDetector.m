function [rate,PF,PD] = MapDetector(threshold,A,sigma,P0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
   P0 = .8;
end
P = P0*100;

iter = 500;
Ri=zeros(iter,1);PFi=Ri;PDi=Ri;

for i=1:iter
% Create Target
target = randi(100,1000,1);
target(target<=P)=0;
target(target>P)=1;

% Add Noise
y = target;
y(target==0) = normrnd(0,sigma,size(y(y==0)));
y(target==1) = normrnd(A,sigma,size(y(y==1)));

% MAP Detection 
d = y;
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