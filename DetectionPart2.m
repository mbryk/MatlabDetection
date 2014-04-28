% Mark Bryk and Yaron Tokayer
% ECE 302
% Detection Exercise Part 2
% May 1, 2014

%% 
clc, clear, close all

%% Part 2 - Photon Detection

P0=.8; P1=1-P0;
lambda0 = 75; lambda1 = 90;
C = [0 1; 1 0];
eta = (C(2,1)-C(1,1))/(C(1,2)-C(2,2)) * (P0/P1);
threshold = (log(eta)+(lambda1-lambda0))/log(lambda1/lambda0);

[rate,PD,PF] = PhotonDetector(threshold,lambda0,lambda1);

ratios = 1:.1:1.5;
thresholds = 70:2:150;
styles = ['b','g','k','c','m','r'];
legends = cell(length(ratios),1);
for i = 1:length(ratios)
    lambda2 = lambda1*ratios(i);
    for j=1:length(thresholds)
        threshold = thresholds(j);
        [rat,PF(i,j),PD(i,j)]=PhotonDetector(threshold,lambda1,lambda2);
    end
    hold on
    plot(PF(i,:),PD(i,:),styles(i));
    hold off
    legends{i} = strcat('\lambda_2/\lambda_1=',num2str(ratios(i)));
end
legend(legends);
title('ROC for Different Photon Rates')
ylabel('P_D, Probability of Detection')
xlabel('P_F, Probability of False Alarm')