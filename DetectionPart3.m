% Mark Bryk and Yaron Tokayer
% ECE 302
% Detection Exercise Part 3
% May 1, 2014

%% 
clc, clear, close all
load('Iris.mat')
% Use first half of each data set for training, and second half for testing.

%% Training Stage
trainData = zeros(25, 4, 3); % preallocate
trainData(:,:,1) = features(1:25, :);
trainData(:,:,2) = features(51:75, :);
trainData(:,:,3) = features(101:125, :);

testInd = [26:50,76:100,126:150];
testData = features(testInd,:);
trueLabels = labels(testInd);

% Create mean and covariance matrices:
mu = zeros(3,4); covMat = zeros(4,4,3); % preallocate
mu(1,:) = mean(trainData(:,:,1)); covMat(:,:,1) = cov(trainData(:,:,1));
mu(2,:) = mean(trainData(:,:,2)); covMat(:,:,2) = cov(trainData(:,:,2));
mu(3,:) = mean(trainData(:,:,3)); covMat(:,:,3) = cov(trainData(:,:,3));

%% Testing Stage

guess = zeros(75,1);

p1 = mvnpdf(testData, mu(1,:), covMat(:,:,1)); % Probabilities of belonging
p2 = mvnpdf(testData, mu(2,:), covMat(:,:,2)); % to a particular class   
p3 = mvnpdf(testData, mu(3,:), covMat(:,:,3));

ps = [p1,p2,p3];
[y,Classes] = max(ps,[],2);
probErr = sum(Classes~=trueLabels)/length(Classes);

confMat = zeros(3);
classInd = 1:25:51; % indeces beginning each class
for i=1:3
    index = classInd(i);
    for j=1:3
        confMat(i,j) = sum(Classes(index:index+24)==j);
    end
end
