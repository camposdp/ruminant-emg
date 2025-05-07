% Daniel Prado de Campos - UTFPR - 03/05/2019
% Processamento de EMG - Parte 2
% Treinamento kSVD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Etapas
% 0) Carregamento dos dados pré-definições
% 1) Separação das instâncias e segmentação
% 2) Definição de kSVD



% *) Funções

%% 0) Carregamento dos dados pré-definições

clear
clc
close all

load EMG_labels_and_windows
load labels_val_e_test


fs = 2000;
N = length(wTv);
Nb = 0.2*fs;
Nf = 0.6*fs;
w0 = Nb+Nf;


d = [1 2 4 8];
k = [0.5 1.0 1.5];



for l = 1:length(d)
for m = 1:length(k)

    
% 1) Separação das instâncias e segmentação

w = w0/d(l);
EMG = [];
EMG = zeros(N,w);
for i = 1:N
EMG(i,:) = decimate(EMGn(xiv(i)-Nb:xiv(i)+Nf-1),d(l)); 
end


% 2) Definição de kSVD

s = size(EMG);
W = s(2);

K=W*k(m);
         
L = K*0.01;

param.L = L;   % number of elements in each linear combination.
param.K = K; % number of dictionary elements
param.numIteration = 10; % number of iteration to execute the K-SVD algorithm.
param.errorFlag = 0; % decompose signals until a certain error is reached. do not use fix number of coefficients.
%param.errorGoal = sigma;
param.preserveDCAtom = 0;
param.InitializationMethod =  'DataElements';
param.displayProgress = 1;


[D{m,l},output{m,l}]  = KSVD(EMG',param);


[A{m,l}]=OMP(D{m,l},EMG',L); 
 

%Y=T*[1 2]';

%mdl_KNNeq{m,n}=fitlda(X,Y);
% e(m,n,:) = kfoldLoss(crossval(mdl_KNNeq),'Mode','individual');



end
end
