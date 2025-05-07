
% Daniel Prado de Campos - UTFPR - 03/05/2019
% Processamento de EMG - Parte 3
% Treinamento DLSI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Etapas
% 0) Carregamento dos dados pr�-defini��e
% 1) Preparar Matriz de sEMG




%% 0) Carregamento dos dados pr�-defini��es

clear
clc
close all


load Segments_labels_Cpartition

fs = 2000;

%label_test = wTt;
%label_train = wTv;

sparsitythres = 5;
d = 1;
lambda = 0.01;
eta = 0.1;

Ki = [100 200 300 400 500 600];

%% 1) Preparar Matriz de EMG

%Usar s� wTv e Y_train (o test ser� usado a posteriori)
for m = 1:length(Ki)
    for l = 1:10

k = Ki(m);
cTrain = Ck.training(l);
label_train = wTv(cTrain);
Y_train_k = Y_train(:,cTrain);

cTest =  Ck.test(l);
label_test = wTv(cTest);
Y_test_k = Y_train(:,cTest);

[acc(l,m), rt] = DLSI_wrapper(Y_train_k, label_train, Y_test_k, label_test, ...
                            k, lambda, eta)
                                                
label_train = [];            
label_test = [];    
Y_train_k = [];  
Y_test_k = [];
cTest = [];
cTrain = [];
c=[];                     

    end
end