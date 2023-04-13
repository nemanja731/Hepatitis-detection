clc, clear, close all

%% Separation of classes
% load haberman.data
T = readtable('hepatitis.csv','TreatAsEmpty',{'?'});
X = T{:, :};
training = X'; 
output = training(end,:);
K1 = training(1:end-1,output == 2);
K2 = training(1:end-1,output == 1);

%% Division into training, validation and test set
N1 = length(K1);
N2 = length(K2);

%% Division into training, validation and test set
K1training = K1(:,1:round(0.6*N1));
K1val = K1(:,round(0.6*N1) + 1 : round(0.8*N1));
K1test = K1(:,round(0.8*N1)+1:N1);

K2training = K2(:,1:0.6*N2);
K2val = K2(:,0.6*N2 + 1: 0.8*N2);
K2test = K2(:,0.8*N2 + 1:N2);

%%
inputTraining = [K1training,K2training];
outputTraining = [ones(1,42),zeros(1,51)];

inputVal = [K1val,K2val];
outputVal = [ones(1,14),zeros(1,17)];

%% 
indTest = randperm(length(outputTraining));
inputTraining = inputTraining(:,indTest);
outputTraining = outputTraining(indTest);

indVal = randperm(length(outputVal));
inputVal = inputVal(:,indVal);
outputVal = outputVal(indVal);

inputAll = [inputTraining,inputVal];
outputAll = [outputTraining,outputVal];
N = length(outputAll);

%% Creation of Neural Network
architecture = [1];
net = patternnet(architecture);

net.trainParam.epochs = 1000;
net.trainParam.goal = 1e-4;
net.trainParam.min_grad = 1e-5;

net.trainFcn = 'trainrp';

net.performFcn = 'mse'; %'crossentropy'; 
 
net.divideFcn = 'divideint';
                
net.divideParam.trainInd = 1:length(inputTraining); %80 percent for training
net.divideParam.valInd =  length(inputTraining) + 1: length(inputTraining) + length(inputVal);%10 percent for validation
net.divideParam.testInd = []; 

net = train(net,inputAll,outputAll);
