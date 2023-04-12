clc, clear, close all

%% Odvajanje klasa
% load haberman.data
T = readtable('hepatitis.csv','TreatAsEmpty',{'?'});
X = T{:, :};
trening = X'; 
izlaz = trening(end,:);
K1 = trening(1:end-1,izlaz == 2);
K2 = trening(1:end-1,izlaz == 1);

%% Podela na trening,validacioni i test skup
N1 = length(K1);
N2 = length(K2);

%% Podela na trening,validacioni i test skup
K1trening = K1(:,1:round(0.6*N1));
K1val = K1(:,round(0.6*N1) + 1 : round(0.8*N1));
K1test = K1(:,round(0.8*N1)+1:N1);

K2trening = K2(:,1:0.6*N2);
K2val = K2(:,0.6*N2 + 1: 0.8*N2);
K2test = K2(:,0.8*N2 + 1:N2);

%%
ulazTrening = [K1trening,K2trening];
izlazTrening = [ones(1,42),zeros(1,51)];

ulazVal = [K1val,K2val];
izlazVal = [ones(1,14),zeros(1,17)];

%% 
indTest = randperm(length(izlazTrening));
ulazTrening = ulazTrening(:,indTest);
izlazTrening = izlazTrening(indTest);

indVal = randperm(length(izlazVal));
ulazVal = ulazVal(:,indVal);
izlazVal = izlazVal(indVal);

ulazSve = [ulazTrening,ulazVal];
izlazSve = [izlazTrening,izlazVal];
N = length(izlazSve);

%% Kreiranje NM
arhitektura = [1];
net = patternnet(arhitektura);

net.trainParam.epochs = 1000;
net.trainParam.goal = 1e-4;
net.trainParam.min_grad = 1e-5;

net.trainFcn = 'trainrp';

net.performFcn = 'mse'; %'crossentropy'; 
 
net.divideFcn = 'divideint';
                
net.divideParam.trainInd = 1:length(ulazTrening); %80 posto za trening
net.divideParam.valInd =  length(ulazTrening) + 1: length(ulazTrening) + length(ulazVal);%10 posto za validaciju
net.divideParam.testInd = []; 

net = train(net,ulazSve,izlazSve);
