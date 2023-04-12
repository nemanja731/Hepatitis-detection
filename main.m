clear
close all
clc

%% Dopuna nedostajucih vrednosti u bazi i prebacivanje csv fajla u matricu X

T = readtable('hepatitis.csv','TreatAsEmpty',{'?'});
[numSamples, numAllAttributes] = size(T);
numAllAttributes = numAllAttributes - 1;
%zamena prve kolone(klasa) sa poslednjom(atribut) radi lakse analize
for j = 1:numSamples
    t = T{j, 1};
    T{j, 1} = T{j, 20};
    T{j, 20} = t;
end
means = zeros(1, numAllAttributes);
for i = 1:numAllAttributes
    numNaNs = 0;
    for j = 1:numSamples
        if isnan(T{j, i})
            numNaNs = numNaNs + 1;
        else
            means(i) = means(i) + T{j, i};
        end
    end
    means(i) = round(means(i)/(numSamples - numNaNs));
end
for i = 1:numAllAttributes
    for j = 1:numSamples
        if isnan(T{j, i})
            if sum(i == [2, 15, 16, 17, 18, 19]) == 0
                T{j, i} = means(i);
            else
                d = rand;
                if d > 0.5
                    k = -1;
                else
                    k = 1;
                end
                T{j, i} = round(means(i) + k*d*means(i)/10);
            end
        end
    end
end
X = T{:, :};

%% Prikaz klasa

numAttributes = 10;
figure
    histogram(X(:, end),'Normalization','probability');
    title('Die                                                Live');
    xlim([0.49 2.51])
    
%% Pearsonov i Spermanov koeficijent korelacije

R1 = corrcoef(X);
[rho, indexesPearson] = sort(R1(end, 1:end - 1),'descend');
% Kreacija promenljive ranga
rank = sort(X);
for i = 1:numAllAttributes + 1
    values = unique(X(:, i)); 
    for rg = 1:length(values)
        rank(X(:, i) == values(rg), i) = rg;
    end
end
R2 = corrcoef(rank);
figure
    heatmap(R1);
    title('Pearsonov koeficijent korelacije');
figure
    heatmap(R2);
    title('Spermanov koeficijent korelacije');
    
%% Informaciona dobit IG

% Prikaz histograma za svaki atribut, za odredjivanje
% kvantizacionih nivoa kod kontinualnih atributa
%{
for i = 1:19
    figure;
        hist(X(:, i))
end
%}
p1 = sum(X(:, end) == 1)/numSamples;
p2 = 1 - p1;
InfoD = -p1*log2(p1) - p2*log2(p2);
InfoDA = zeros(1, numAllAttributes);
type = ['D  ';'AGE';'D  ';'D  ';'D  ';'D  ';'D  ';'D  ';'D  ';'D  ';'D  ';'D  ';'D  ';'D  ';'BIL';'ALK';'SGO';'ALB';'PRO'];
for i = 1:numAllAttributes
    InfoDA(i) = infoDAFunction(InfoD, X(:, i), X(:, end), type(i, :));
end
[~, indexesInfoDA] = sort(InfoDA, 'descend');
Z = X(:, sort(indexesInfoDA(1:numAttributes)));

%% Podela na obucavajuci i testirajuci skup

X1 = Z(X(:, end) == 1, :);
X2 = Z(X(:, end) == 2, :);
indX1 = randperm(length(X1));
indX2 = randperm(length(X2));
percentage = 0.7;
% Obucavajuci skup
X1o = X1(indX1(1:floor(percentage*end)), :);
X2o = X2(indX2(1:floor(percentage*end)), :);
% Testirajuci skup
X1t = X1(indX1(floor(percentage*end) + 1:end), :);
X2t = X2(indX2(floor(percentage*end) + 1:end), :);

%% Redukcija dimenzija na bazi matrica rasejanja

M1 = mean(X1o)';
M2 = mean(X2o)';
M0 = M1*p1 + M2*p2;
S1 = cov(X1o);
S2 = cov(X2o);
% Matrice
Sw = S1*p1 + S2*p2;
Sb = (M1 - M0)*(M1 - M0)'*p1 + (M2 - M0)*(M2 - M0)'*p2;
Sm = Sb + Sw;
T = Sw^-1*Sb;

%% 2 dimenzije

[V, ~] = eigs(T, 2);
A = real(V);
Y1o = A'*X1o';
Y2o = A'*X2o';
Y1t = A'*X1t';
Y2t = A'*X2t';
figure;
    plot(Y1o(1, :), Y1o(2, :), 'ko', Y2o(1, :), Y2o(2, :), 'rx');
    title('LDA obucavajuceg skupa na 2 dimenzije')
figure;
    plot(Y1t(1, :), Y1t(2, :), 'ko', Y2t(1, :), Y2t(2, :), 'rx');
    title('LDA testirajuceg skupa na 2 dimenzije')
    
%% 3 dimenzije

[V, ~] = eigs(T, 3);
A = real(V);
Y1o = A'*X1o';
Y2o = A'*X2o';
Y1t = A'*X1t';
Y2t = A'*X2t';
figure;
    plot3(Y1o(1, :), Y1o(2, :), Y1o(3, :), 'ko', Y2o(1, :), Y2o(2, :), Y2o(3, :), 'rx');
    title('LDA obucavajuceg skupa na 3 dimenzije')
    grid on;
figure;
    plot3(Y1t(1, :), Y1t(2, :), Y1t(3, :), 'ko', Y2t(1, :), Y2t(2, :), Y2t(3, :), 'rx');
    title('LDA testirajuceg skupa na 3 dimenzije')
    grid on;
    
%% Testiranje hipoteza

% skupovi su nam nebalansirani
c = zeros(2,2);
a = sum(X(:,end) == 1);
b = sum(X(:,end) == 2);
c11 = 0;
c22 = 0;
c21 = 1;
c12 = 10;
% za sad najbolje 1,10
prag = -log((c22-c12)/(c11-c21));

for i=1:length(X1t)
 % prolazimo kroz klasu X1_t i odredjujemo h(X1)
 f1=(1/(2*pi*(det(S1))^0.5))*exp(-0.5*(X1t(i,:)'-M1)'*inv(S1)*(X1t(i,:)'-M1));
 f2=(1/(2*pi*(det(S2))^0.5))*exp(-0.5*(X1t(i,:)'-M2)'*inv(S2)*(X1t(i,:)'-M2));
 h=-log(f1)+log(f2);
 if h<prag % klasifikator donosi odluku da je klasa 1
    c(1,1)=c(1,1)+1;
 else
    % klasifikator gresi i misli da je klasa 2
    c(2,1)=c(2,1)+1;
 end
end

for i=1:length(X2t)
 % prolazimo kroz klasu X2t i odredjujemo h(X2)
 f1=(1/(2*pi*(det(S1))^0.5))*exp(-0.5*(X2t(i,:)'-M1)'*inv(S1)*(X2t(i,:)'-M1));
 f2=(1/(2*pi*(det(S2))^0.5))*exp(-0.5*(X2t(i,:)'-M2)'*inv(S2)*(X2t(i,:)'-M2));
 h=-log(f1)+log(f2);
 if h<prag % klasifikator gresi i misli da je klasa 1
    c(1,2)=c(1,2)+1;
 else
    % klasifikator doneo ispravnu odluku
    c(2,2)=c(2,2)+1;
 end
end

err = (c(1,2)+c(2,1))/sum(sum(c)) 
c

podaci = X;

%% Transponovanje matrice i iscrtavanje klasa
ulaz = podaci(:, 1:19)';
klasa = podaci(:, 20)';
N = length(klasa);

K1 = ulaz(:, klasa == 1);
K2 = ulaz(:, klasa == 2);

izlaz = zeros(2, N);

izlaz(1, klasa == 1) = 1;
izlaz(2, klasa == 2) = 1;

%% Podela na trening i test skup  
ind = randperm(N);
indTrening = ind(1 : round(0.7*N));
indTest = ind(round(0.7*N+1) : N);

ulazTrening = ulaz(:, indTrening);
izlazTrening = izlaz(:, indTrening);

ulazTest = ulaz(:, indTest);
izlazTest = izlaz(:, indTest);

%% Kreiranje NM
arhitektura = [10 10 10];
net = patternnet(arhitektura);

net.trainParam.epochs = 1000; %max epoha
net.trainParam.goal = 1e-4; %ciljana greska
net.trainParam.min_grad = 1e-5;

net.trainFcn = 'trainbfg';

net.divideFcn = ''; 
net.performFcn = 'msereg'; %'crossentropy'; 
net.performFcn = 'crossentropy'; %'crossentropy'; 
net.performParam.regularization = 5e-5;
%% Treniranje NM
net = train(net, ulazTrening, izlazTrening);
%% Performanse NM
predTest = net(ulazTest);
predTrening = net(ulazTrening);
figure, plotconfusion(izlazTest, predTest)
figure, plotconfusion(izlazTrening, predTrening)