%% MATLAB Detection Exercise
%  Mark Bryk and Yaron Tokayer
%  ECE 302 - Stochastics and Probability
%  5/1/14
%

%% 
clc, clear, close all

%% Part 1 - Radar Detection

% A 
C = [0 1; 1 0]; P0 = .8; A = 5; sigma = 2;

P1 = 1-P0;
eta = (C(2,1)-C(1,1))/(C(1,2)-C(2,2)) * (P0/P1);
gamma = A/2 + (sigma^2)*log(eta)/A;
[rateA,PF,PD] = MapDetector(gamma,A,sigma);

%% B
snr = [.1,.25,1,2,4]; 
thresholds = -5:.2:10;
styles = ['b','g','k','y','m'];
sigmas = sqrt(A./snr);
legends = cell(length(sigmas),1);
for i = 1:length(sigmas)
    sigma = sigmas(i);
    for j=1:length(thresholds)
        threshold = thresholds(j);
        [rat(i,j),PF(i,j),PD(i,j)]=MapDetector(threshold,A,sigma);
    end
    hold on
    plot(PF(i,:),PD(i,:),styles(i));
    hold off
    legends{i} = strcat('SNR=',num2str(snr(i)));
end
legend(legends);

%% C - For SNR = 1
C = [0 10; 1 0]; P0 = .8; A = 5; 
index = 3;
snrC = snr(index); sigma = sqrt(A/snrC);
P1 = 1-P0;
eta = (C(2,1)-C(1,1))/(C(1,2)-C(2,2)) * (P0/P1);
gamma = A/2 + (sigma^2)*log(eta)/A;
[rateC,PFC,PDC] = MapDetector(gamma,A,sigma);
% [t_delta, t_ind] = min(abs( PF(index,:)-PFC ));
hold on
% plot(PF(index,t_ind),PD(index,t_ind),'r*','MarkerSize',8);
plot(PFC,PDC,'r*','MarkerSize',8);
hold off

%% D - For SNR = 1
P1s = .1:.02:.9;
C = [0 1; 1 0]; A = 5; 
snr = 1; sigma = sqrt(A/snr);
cost = zeros(length(P1s),1);
for i=1:length(P1s)
    P1 = P1s(i);
    P0 = 1-P1;
    eta = (C(2,1)-C(1,1))/(C(1,2)-C(2,2)) * (P0/P1);
    gamma = A/2 + (sigma^2)*log(eta)/A;
    [rateD,PF,PD] = MapDetector(gamma,A,sigma,P0);
    PM = 1-PD;
    cost(i) = PF*C(2,1)*P0 + PM*C(1,2)*P1;
end
figure,plot(P1s,cost);
