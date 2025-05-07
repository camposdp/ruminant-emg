% Daniel Prado de Campos - UTFPR - 03/05/2019
% Processamento de EMG - Parte 3
% Treinamento LC-kSVD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Etapas
% 0) Carregamento dos dados pré-definiçõe
% 1) Preparar Matriz de sEMG
% 2) Parâmetros LC-kSVD
% 3) Aprendizado do dicionário 
% 4) LC-kSVD1
% 5) LC-kSVD2



%% 0) Carregamento dos dados pré-definições


clear
clc
close all

%load Dic_kSVD
load EMG_labels_and_windows
load labels_val_e_test


WTv = wTv;
WTt = wTt;

Xiv = xiv;
Xit = xit;

fs = 2000;

Nb = 0.2*fs;
Nf = 0.6*fs;
w0 = Nb+Nf;




Nt = length(wTt);
Nv = length(wTv);

V = [500 1000 1500 2000 3000 4000];

for l = 1:length(V)
for n = 1:2

    
    
%% 1) Preparar Matriz de EMG

d = 8;
w = w0/d;


N_v = V(l);
N_t = V(l)*0.2;

v=sort(randperm(Nv,N_v));
t=sort(randperm(Nt,N_t));

Y_train = zeros(w,N_v);
Y_test = zeros(w,N_t);


wTv = WTv(v);
wTt = WTt(t);

xiv = Xiv(v);
xit = Xit(t);



label_train = wTv;
for i = 1:N_v
Y_train(:,i) = decimate(EMGn(xiv(i)-Nb:xiv(i)+Nf-1),d); 
end

if n == 1
    label_test = wTt;
    
    for i = 1:N_t
    Y_test(:,i) = decimate(EMGn(xit(i)-Nb:xit(i)+Nf-1),d); 
    end

elseif n == 2
    Y_test = Y_train;
    label_test = label_train;
end

%% 3) Parâmetros LC-kSVD

%ki = [500 1000 2000];
%spi = [5 10 15 20 30 50];



%range_train = label_to_range(label_train);
%range_test = label_to_range(label_test);


k = 50;
sparsitythres  = 5;

valpha = 0.01;
vbeta = 0.01;

 acc = zeros(1, 2);
rt = zeros(1, 2);
    C = max(label_train);
    H_train = lcksvd_buildH(label_train);
    H_test = lcksvd_buildH(label_test);


    sqrt_alpha = sqrt(valpha); % weights for label constraint term
    sqrt_beta = sqrt(vbeta); % weights for classification err term
    dictsize = C*k; % dictionary size
    iterations = 50; % iteration number
    iterations4ini = 30; % iteration number for initialization

%% 4) Aprendizado do dicionário 
    % get initial dictionary Dinit and Winit
    fprintf('\nLC-KSVD initialization... ');
    [Dinit,Tinit,Winit,Q_train] = initialization4LCKSVD(Y_train, H_train, ...
        dictsize, iterations4ini, sparsitythres);
    fprintf('done!');
    
% %% 5) ========= LCKSVD1 ==============================  
%     % run LC K-SVD Training (reconstruction err + class penalty)
    fprintf('\nDictionary and classifier learning by LC-KSVD1...');
    tic;
    [D1,X1,T1,W1] = labelconsistentksvd1(Y_train, Dinit, Q_train, Tinit, ...
        H_train,iterations,sparsitythres,sqrt_alpha);
    rt(1) = toc;
    fprintf('done!');
%     % classification process
%     [prediction1, acc(1)] = classification(D1, W1, Y_test, H_test, sparsitythres);
     %fprintf('\nFinal recognition rate for LC-KSVD1 is : %.03f \n', acc(1));
%% 6) ========= LCKSVD2 ==============================  
    % run LC k-svd training (reconstruction err + class penalty + classifier err)
%     fprintf('\nDictionary and classifier learning by LC-KSVD2...');
%     tic;
%     [D2,X2,T2,W2] = labelconsistentksvd2(Y_train, Dinit, Q_train, ...
%         Tinit, H_train, Winit, iterations, sparsitythres, sqrt_alpha, sqrt_beta);
%     rt(2) = toc;
%     fprintf('done!\n');
    %classification
  %  [prediction2, acc(2)] = classification(D2, W2, Y_test, H_test, sparsitythres);
%     
    X = myOMP(Y_test, D1, sparsitythres);
    Mdl = fitcdiscr(full(X)',label_test,'DiscrimType','pseudolinear');
    e = kfoldLoss(crossval(Mdl));
    acc(1) = 1-e; 
 %   fprintf('\nFinal recognition rate for LC-KSVD2 is : %.03f \n', acc(2));
 
  fprintf('\nFinal recognition rate for LC-KSVD1 is : %.03f \n', acc(1));
    

     ACC(l,n) = acc(1);
%    ACC2(l,n) = acc(2);
    
    
end
end

% subplot(1,2,1)
% plot(V,ACC)
% ylabel('Acc')
% legend('test','train')
% subplot(1,2,2)
plot(V,ACC2)
ylabel('Acc')
