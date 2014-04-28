function [ rate,PF,PD ] = PhotonDetector(threshold,lambda0,lambda1)

iter = 100;
Ri=zeros(iter,1);PFi=Ri;PDi=Ri;
for i=1:iter
    
% Create Signal
bits = randi(10,1000,1);
bits(bits<=8)=0;
bits(bits>8)=1;

% Send Photons
y = bits;
y(bits==0) = poissrnd(lambda0,size(y(bits==0)));
y(bits==1) = poissrnd(lambda1,size(y(bits==1)));

% LRT Comparison
d=y;
d(y<=threshold) = 0;
d(y>threshold) = 1;

% Calculate Statistics
Ri(i) = sum(bits~=d)/length(bits);
H0 = d(bits==0);
H1 = d(bits==1);
PFi(i) = sum(H0)/length(H0);
PDi(i) = sum(H1)/length(H1);
end

rate = mean(Ri);
PF = mean(PFi);
PD = mean(PDi);
end