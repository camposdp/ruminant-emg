% Daniel Prado de Campos - UTFPR - 03/05/2019
% Processamento de EMG - Parte 3
% LC-kSVD - Sparse Discriminant Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Etapas
% 0) Carregamento dos dados pré-definiçõe
% 1) Preparar Matriz de sEMG
% 2) LC-kSVD
% 3) SDA



%% 0) Carregamento dos dados e pré-definições
% 

close all
clear 
clc

% % % load EMG_raw_labels
% % % load EMGre_exp
% % % load Segments_labels_Cpartition
% % % load Segments_labels_index
% % % 
% % % fsa = 2000; %frequência de amostragem
% % % 
% % % %% 1) Preparar Matriz de EMG
% % % 
% % % BL = [4.45 4.7]*1e5;
% % % 
% % % Wcrit = [0.1*fsa 1*fsa];
% % % Wsmooth = 0.15*fsa;
% % % k=4.0;
% % % 
% % % x = smoothdata(abs(EMGn'),'movmean',Wsmooth);
% % % th = mean(EMGt_re(BL(1):BL(2)))+k*std(EMGt_re(BL(1):BL(2)));
% % % 
% % % [xi,xf,wT] = DTSeg(EMGn',Wcrit,x,th,T);
% % % 
% % % xi = xi(1:end-1);
% % % xf = xf(1:end-1);
% % % wT = wT(1:end-1);
% % % 
% % % %Extrair features
% % % idx = find(ind_train==1);
% % % %Usar só wTv e Y_train (o test será usado a posteriori)
% % % wTx = wT(idx);
% % % xix = xi(idx);
% % % YT = wTx;
% % % 
% % % lambda1 = 0.01;
% % % lambda2 = 0.01;
% % % 
% % % 
% % % % Wa = [0 0.05 0.1 0.2 0.3];
% % % % Wb = [0.3 0.35 0.4 0.45 0.5 0.6];
% % % %best: W = [0 0.5]
% % % % for i = 1:length(Wa)
% % % % for j = 1:length(Wb)
% % % % 
% % % % W = [Wa(i)*fsa Wb(j)*fsa];
% % % 
% % % W = [0.2 0.4]*fsa;
% % % 
% % % Y_n = [];
% % % for ind = 1:length(idx)
% % % Y_n(:,ind) = EMGn(xix(ind)-W(1):xix(ind)+W(2)-1);
% % % end
% % % 
% % % %% 2) LC-kSVD
% % % 

% % %         
% % % i=1;        
% % % cTrain = Ck.training(i);
% % % label_train = wTx(cTrain);
% % % Y_train = Y_n(:,cTrain);
% % % 
% % % cTest =  Ck.test(i);
% % % label_test = wTx(cTest);    
% % % Y_test = Y_n(:,cTest);


%%%%%%%%%%%%%%%%%%%%%%%%%%

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%


k = 50;
sparsitythres = 20;
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

    
    Xtr = myOMP(Y_train, D1, sparsitythres);
    Xtt = myOMP(Y_test, D1, sparsitythres);
    Ytr = label_train';
    Ytt = label_test';
    %% SDA
    
    
    % set parameters
lambda = 1e-6; % l2-norm
stop = [-1]; % l1-norm. negative: number of vars in LARS-EN
maxiter = 25; % max iter in SDCA alg.

I = eye(3);
Ytr_m = I(Ytr,:);



% perform SDA
[sl theta] = slda(full(Xtr)', Ytr_m, lambda, stop, 10,maxiter, 1e-6,1);
% Project data onto the sparse directions (dim=2)
DCtr = full(Xtr)'*sl;
DCtt = full(Xtt)'*sl;

% Classification (LDA of projected data)
[class,trerr] = classify(DCtr,DCtr,Ytr,'linear');
trerr = length(find(class~=Ytr))/length(Ytr)

[class,tsterr] = classify(DCtt,DCtr,Ytr,'linear');
tsterr = length(find(class~=Ytt))/length(Ytt)

% plot sparse discriminative directions
% figure;
% plot(DC(I0,1),DC(I0,2),'ro','linewidth',2,'MarkerSize',9), hold on
% plot(DC(I1,1),DC(I1,2),'ks','linewidth',2,'MarkerSize',9)
% plot(DC(I2,1),DC(I2,2),'bv','linewidth',2,'MarkerSize',9)
% xlabel('1st SD'), ylabel('2nd SD')
% legend('Mel','Pol','Ven',4)

    






function [xi,xf,wT] = DTSeg(data,Wcrit,x,th,T)
%Double Threshold segmentation (variable length)
% A.K.A Onset Segmentation
%
% -> Find segments which are:          
%    *Longer than Wcrit               
%    *Greater than th                  
%
% Inputs:
% - data: data vector
% - win: Signal smooth window size
% - Wcrit: Time threshold (min segmentation length)
% - th: Amplitude threshold
% - method: smoothing method
% Outputs:
% - segM: segment cell
% - where: segment extension

%data length
N= max(size(data));
%smooth absolute signal

%PB = fH*2/fs;
%[Bb,Ab] = butter(6,PB,'low');
%x=filter(Bb,Ab,abs(ex));


%th = mean(x(BL(1):BL(2)))+k*std(x(BL(1):BL(2)));

%th = k*std(x(BL(1):BL(2)));

%segM=[];
%where=zeros(1,N);
wT = [];

i=1;
j=1;
while i~=N

   %Signal is under th... 
   while abs(x(i)<th) && i~=N 
       i=i+1;        
   end
   %...signal is over th!
   
   %Since Ci, there is a new segment...
   Ci=i;
   
   %...still over th....
   while abs(x(i)>th) && i~=N  
       i=i+1;       
   end
   %...and its gone. Now it is under th.
    
   %The segment ends at Cf. 
   Cf=i;
   
   %The segment existance is limited by Ci and Cf.
   %The difference is the segment length 
   C=Cf-Ci;
   
   %Is this segment longer than Wcrit?
   if (C>Wcrit(1)) && (C<Wcrit(2))
       
       xi(j) = Ci;
       xf(j) = Cf;
       %We got a new segment!
       w=Ci:Cf;             %segment boundary
       %w_fix=Ci:(Ci+W-1)
       
       Tw = T(w,:);
       [~,wi]=max(sum(Tw));
       wT=cat(2,wT,wi);
       
       %where(w) = 1;        %If there is a segment, there is a 1.
       %segM(j,:) = data(w);   %Add the segment to the cell
       %segM = data(w);   %Add the segment to the cell
       
       
       
       j=j+1;
       
       
    
   end
   
end

end


    
    
    
    