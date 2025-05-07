



%% 0) Carregamento dos dados pré-definições

clear
clc
close all


load Segments_labels_Cpartition

fs = 2000;

%label_test = wTt;
%label_train = wTv;

sparsitythres = 5;
d = 1;

lambda1 = 0.01;
lambda2 = 0.01;

Ki = [100 200 300 400 500 600 700];

%% 1) Preparar Matriz de EMG

%Usar só wTv e Y_train (o test será usado a posteriori)
for m = 1:length(Ki)
    for l = 1:10

k = Ki(m);
cTrain = Ck.training(l);
label_train = wTv(cTrain);
Y_train_k = Y_train(:,cTrain);

cTest =  Ck.test(l);
label_test = wTv(cTest);
Y_test_k = Y_train(:,cTest);

[AC, rt] = FDDL_wrapper(Y_train_k, label_train, Y_test_k, label_test, ...
                            k, lambda1, lambda2)
           
acc(l,m) = max(AC)
 
label_train = [];            
label_test = [];    
Y_train_k = [];  
Y_test_k = [];
cTest = [];
cTrain = [];
c=[];                     

    end
end