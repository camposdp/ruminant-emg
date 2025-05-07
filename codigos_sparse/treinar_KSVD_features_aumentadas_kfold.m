% Daniel Prado de Campos - UTFPR - 20/02/2019
% Processamento de EMG do bovino para classificação entre
% ruminação/ócio/alimentação
% Comparação com níveis de ruído entre melhores
% resultados de KSVD e handcrafted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0) Carregamento dos dados e pré-definições
% 

close all
clear 
clc

load EMG_raw_labels
load EMGre_exp
load Segments_labels_Cpartition
load Segments_labels_index

fsa = 2000; %frequência de amostragem

%% 3) Segmentação DTOS e extração de características


K = 300;
L = 5;

    param.L = L;   % number of elements in each linear combination.
 % number of dictionary elements
param.numIteration = 30; % number of iteration to execute the K-SVD algorithm.
param.errorFlag = 0; % decompose signals until a certain error is reached. do not use fix number of coefficients.
%param.errorGoal = sigma;
param.preserveDCAtom = 0;
param.InitializationMethod =  'DataElements';
param.displayProgress = 1;
    





BL = [4.45 4.7]*1e5;

Wcrit = [0.1*fsa 1*fsa];
Wsmooth = 0.15*fsa;
k=4.0;

x = smoothdata(abs(EMGn'),'movmean',Wsmooth);
th = mean(EMGt_re(BL(1):BL(2)))+k*std(EMGt_re(BL(1):BL(2)));

[xi,xf,wT] = DTSeg(EMGn',Wcrit,x,th,T);

xi = xi(1:end-1);
xf = xf(1:end-1);
wT = wT(1:end-1);

W = [0.0 0.6]*fsa;




%% 1) Preparar Matriz de EMG


%Extrair features
idx = find(ind_train==1);
%Usar só wTv e Y_train (o test será usado a posteriori)
wTx = wT(idx);
YT = wTx;


param.L = 2;


%snr = [-10 0 5 10 20];

Ki = [300];
for i = 1:length(Ki)
    K = Ki(i);
param.K = K;

    for l = 1:10

        
% 
%         for s = 1:2400
%         Y_n(:,s) = awgn(Y_train(:,s),snr(i),'measured');
%         end
Y_n = Y_train;
       
        
        
        cTrain = Ck.training(l);
        label_train = wTv(cTrain);
        Y_train_k = Y_n(:,cTrain);
        
        cTest =  Ck.test(l);
        label_test = wTv(cTest);    
        Y_test_k = Y_n(:,cTest);
        
        [Di,~] = KSVD(Y_train_k,param);
         D = Di(:,randperm(K));
         [A]=OMP(D,Y_train_k,L);
         %F_train = [sum(abs(A))' max(abs(A))'];
         F_train=aumFeat(A,D);

         mdl_LDA=fitcdiscr(F_train',label_train);
        

        [As2]=OMP(D,Y_test_k,L);
         A2 = full(As2);
         F_test=aumFeat(A2,D);
         %F_test = [sum(abs(A2))' max(abs(A2))'];
         e(i,l)=loss(mdl_LDA,F_test',label_test)

        
        
        

    end


end





%% Funções
% 




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
% 
% 
% function F = aumFeat(A,D)
% 
% [~,idx]=sort(abs(A(find(A(:,2)),2)));
% L=length(idx);
% N=length(A);
% F=[];
% for k = 1:N
% x=(A(find(A(:,k)),k).*D(:,find(A(:,k)))')';
% 
% for j = 1:L
% i = idx(j);
% 
% X = full(x);
% Sf = X(:,i);
% 
% th_ssc=0.001;
% th_zc=0.11;
% th_wamp = 0.06;
% 
% Ff(1)=ZC(th_zc,Sf);
% Ff(2)=WAMP(th_wamp,Sf);
% Ff(3)=RMS(Sf);
% Ff(4)=VAR(Sf);
% Ff(5)=DASDV(Sf);
% Ff(6)=IAV(Sf);
% Ff(7)=MFL(Sf);
% Ff(8)=MSR(Sf);
% Ff(9)=LS(Sf,2);
% 
% F = [F;Ff];
% end
% Fi(:,k) = F;
% F = [];
% end
% 
% end


 function F=pooling(A)
F = [];
 K = [1 2 4 5 10];
 for j = 1:5
 k = K(j);
% A = full(As);
% k=4;
S = size(A);
N = S(1);
%P=permn([1 0],k);
%P = P(1:end-1,:);
P = eye(k);
n = N/k;

for i = 1:k
    ind = [1+n*(i-1):n*i];
    f1(i,:) = sum(abs(A(ind,:)));
    f2(i,:) = max(abs(A(ind,:)));
end

 for q = 1:length(P)
     F1(q,:) = sum(f1(find(P(q,:)),:),1);
     F2(q,:) = max(f2(find(P(q,:)),:),[],1);
 end

F = [[F1' F2'] F];

%F = [f1' f2'];
end
end