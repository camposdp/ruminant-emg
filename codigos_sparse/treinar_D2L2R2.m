
% Daniel Prado de Campos - UTFPR - 03/05/2019
% Processamento de EMG - Parte 3
% Treinamento D2L2R2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Etapas
% 0) Carregamento dos dados pré-definiçõe
% 1) Preparar Matriz de sEMG




%% 0) Carregamento dos dados pré-definições

clear
clc
close all

%load Dic_kSVD
load EMG_labels_and_windows
load labels_val_e_test


fs = 2000;

Nb = 0.2*fs;
Nf = 0.6*fs;
w0 = Nb+Nf;

label_test = wTt;
label_train = wTv;


N_t = length(label_test);
N_v = length(label_train);


%% 1) Preparar Matriz de EMG

d = 1;
w = w0/d;
Y_train = zeros(w,N_v);
Y_test = zeros(w,N_t);

for i = 1:N_v
Y_train(:,i) = decimate(EMGn(xiv(i)-Nb:xiv(i)+Nf-1),d); 
end

for i = 1:N_t
Y_test(:,i) = decimate(EMGn(xit(i)-Nb:xit(i)+Nf-1),d); 
end


lambda1 = 0.01;
lambda2 = 0.01;
alpha = 0.01;
k = 10;




[acc, rt] = D2L2R2_wrapper(Y_train, label_train, Y_test, label_test,...
                        k, lambda1, lambda2, alpha)