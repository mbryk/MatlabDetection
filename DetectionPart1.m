% Mark Bryk and Yaron Tokayer
% ECE 302
% Detection Exercise Part 1
% May 1, 2014

%% 
clc, clear, close all

%% A
C = [0 1; 1 0]; P0 = .8; % MAP cost and given prior
A = 5; sigma = 2; % Choose SNR and mean

P1 = 1-P0;
eta = (C(2,1)-C(1,1))/(C(1,2)-C(2,2)) * (P0/P1);
gamma = A/2 + (sigma^2)*log(eta)/A;
[rateA,PF,PD] = RadarDetector(gamma,A,sigma); % returns rate of error and an estimated prior

%% B
snr = [.1,.25,1,2,4]; 
thresholds = -5:.2:10;
styles = ['b','k','g','m','r'];
sigmas = sqrt(A./snr);
legends = cell(length(sigmas),1);
PF = zeros(length(snr), length(thresholds)); % preallocate memory
PD = zeros(length(snr), length(thresholds));
% Create ROC for each sigma
for i = 1:length(sigmas)
    sigma = sigmas(i);
    for j=1:length(thresholds)
        threshold = thresholds(j);
        [rat(i,j),PF(i,j),PD(i,j)]=RadarDetector(threshold,A,sigma);
    end
    hold on
    plot(PF(i,:),PD(i,:),styles(i));
    hold off
    legends{i} = strcat('SNR=',num2str(snr(i)));
end
title('ROC for Various SNRs')
ylabel('P_D, Probability of Detection')
xlabel('P_F, Probability of False Alarm')
legend(legends);

%% C - For SNR = 1
C = [0 10; 1 0]; P0 = .8; A = 5; % New cost structure
index = 3; % SNR of 1

figROC = figure; plot(PF(3,:), PD(3,:))
title('ROC for SNR = 1')
ylabel('P_D, Probability of Detection')
xlabel('P_F, Probability of False Alarm')

snrC = snr(index); sigma = sqrt(A/snrC);
P1 = 1-P0;
eta = (C(2,1)-C(1,1))/(C(1,2)-C(2,2)) * (P0/P1);
gamma = A/2 + (sigma^2)*log(eta)/A;
[rateC,PFC,PDC] = RadarDetector(gamma,A,sigma);

hold on
plot(PFC,PDC,'r*','MarkerSize',8);
hold off

%% D - For SNR = 1
P1s = 0:.02:1;
C = [0 10; 1 0]; A = 5; 
snr = 1; sigma = sqrt(A/snr);
cost = zeros(length(P1s),1);
for i=1:length(P1s)
    P1 = P1s(i);
    P0 = 1-P1;
    eta = (C(2,1)-C(1,1))/(C(1,2)-C(2,2)) * (P0/P1);
    gamma = A/2 + (sigma^2)*log(eta)/A;
    g(i)=gamma;
    [rateD,PFi,PDi] = RadarDetector(gamma,A,sigma,P0);
    PM = 1-PDi;
    cost(i) = PFi*C(2,1)*P0 + PM*C(1,2)*P1;
end
figure,plot(P1s,cost);
title('E[C](P_1) for C_{01}=10C_{10}')
ylabel('E[C], expecation of cost structure')
xlabel('P_1, probability of a transmission')

%% E
C = [0 10; 1 0];
% Determine the line Pd = a - m*Pf
a = (C(1,2) - C(1,1))/(C(1,2) - C(2,2));
m = (C(2,1) - C(1,1))/(C(1,2) - C(2,2));
PDmm = a - m*PF(3,:);

% Find the intersection point to the ROC by finding the minimum error 
% between elements of PD and PDmm
err = abs(PD(3,:) - PDmm);
[minMSE,minInd] = min(err);
gammaMm = thresholds(minInd);
etaMm = exp(gammaMm*snr);
PFmm = PF(3,minInd); PDmm = PD(3,minInd);
[maxC,estCind] = max(cost);
estC = maxC*ones(1,length(P1s));
P1s(estCind) % minimax estimated P1

% overlay estimate onto graph from part d
figure(gcf), hold on
plot(P1s, estC, 'g')
legend('Known Prior', 'minimax Cost')
hold off

% overlay line onto ROCs
figure(figROC), hold on
plot(PF(3,:), PDmm, 'g')
plot(PFmm,PDmm,'k*','MarkerSize',8);
legend('ROC', '\Gamma for C_{01} = 10C_{10}', 'Minimax P_D line',...
    'location of \Gamma_{minimax}', 'Location', 'SouthEast')
hold off

%% F
C = [0 1; 1 0]; A = 5; sigma = 1;
P0 = .8; P1 = 1-P0;
eta = (C(2,1)-C(1,1))/(C(1,2)-C(2,2)) * (P0/P1);

sigma2 = sigma*10;
ss = sigma^2; s2s = sigma2^2;
gamma = 2*((ss*s2s)/(ss-s2s))*log(eta*sigma2/sigma);
[rateF,PFg,PDg] = RadarDetector2(gamma,A,sigma,sigma2);
    
ratios = [1,2.5,5,10,15];
thresholds = -10:.5:20;
figure, hold on
styles = ['b','k','g','m','r'];
for i=1:length(ratios)
    sigma2 = sigma*ratios(i);
    for j=1:length(thresholds)
        threshold = thresholds(j);
        [rat,PF(j),PD(j)] = RadarDetector2(threshold,A,sigma,sigma2);
    end
    plot(PF,PD,styles(i));
    legends{i} = strcat('\sigma_2/\sigma=',num2str(ratios(i)));
end
title('ROC for Various Variances')
ylabel('P_D, Probability of Detection')
xlabel('P_F, Probability of False Alarm')
legend(legends,'Location','SouthEast');
hold off