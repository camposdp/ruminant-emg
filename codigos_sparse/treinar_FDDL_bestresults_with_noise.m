% Daniel Prado de Campos - UTFPR - 20/02/2019
% Processamento de EMG do bovino para classificação entre
% ruminação/ócio/alimentação
% Comparação com níveis de ruído entre melhores
% resultados de FDDL e handcrafted
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


snr = [-5 0 5 20];



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

%W = [0.05 0.45]*fsa;




%% 1) Preparar Matriz de EMG


%Extrair features
idx = find(ind_train==1);
%Usar só wTv e Y_train (o test será usado a posteriori)
wTx = wT(idx);
xix = xi(idx);
YT = wTx;


lambda1 = 0.01;
lambda2 = 0.01;



% Wa = [0 0.05 0.1 0.2 0.3];
% Wb = [0.3 0.35 0.4 0.45 0.5 0.6];
%best: W = [0 0.5]
% for i = 1:length(Wa)
% for j = 1:length(Wb)
% 
% W = [Wa(i)*fsa Wb(j)*fsa];

W = [0 0.5]*fsa;

Y_train = [];
for ind = 1:length(idx)
Y_train(:,ind) = EMGn(xix(ind)-W(1):xix(ind)+W(2));
end

Ki = [50 100 150 200 250 300 350 400 450 500 550 600];
for j = 1:length(Ki)
Kdic = Ki(j);

%     for j = 1:length(idx)
% 
%     id = idx(j);    
% 
%     Sf = awgn(EMGn(xi(id)-W(1):xi(id)+W(2)),snr(i),'measured');
% 
%     Ff(1,j)=MAV(Sf);
%     Ff(2,j)=WL(Sf);
%     Ff(3,j)=ZC(th_zc,Sf);
%     Ff(4,j)=SSC(th_ssc,Sf);
%     Ff(5,j)=WAMP(th_wamp,Sf);
%     Ff(6,j)=RMS(Sf);
%     Ff(7,j)=VAR(Sf);
%     Ff(8,j)=DASDV(Sf);
%     Ff(9,j)=IAV(Sf);
%     Ff(10,j)=MFL(Sf);
%     Ff(11,j)=MSR(Sf);
%     Ff(12,j)=LS(Sf,2);
% 
% 
% 
%     end
% 
% 
% 
% 
% 
% Xtrainf = Ff;
% 
% 
% 
% %Testar com conjunto de Hudgins
% %YT= wT;
% %Xv = Xtrainv(SET{3},:)';
% Xf = Xtrainf(SET,:)';
% %Tmat = I(YT,:);
% 
% 
% mdl_LDA_Tf = fitcdiscr(Xf,YT);
% cvmodelf = crossval(mdl_LDA_Tf,'CVPartition',Ck);
% enoisy(:,i) = kfoldLoss(cvmodelf,'Mode','individual');
% 
% Xv = [];
% Xf = [];



    for i = 1:10



%         for s = 1:2400
%         Y_n(:,s) = awgn(Y_train(:,s),snr(i),'measured');
%         end
        Y_n = Y_train;
       
        
        
        cTrain = Ck.training(i);
        label_train = wTx(cTrain);
        Y_train_k = Y_n(:,cTrain);
        
        cTest =  Ck.test(i);
        label_test = wTx(cTest);    
        Y_test_k = Y_n(:,cTest);
        
        [AC, rt, Di, CoefMi, optsi] = FDDL_wrapper(Y_train_k, label_train, Y_test_k, label_test, ...
                                    Kdic, lambda1, lambda2);

        acc(i,j) = max(AC)
        D{i,j} = Di;
        CoefM{i,j} = CoefMi;
        opts{i,j}=optsi;
        
        label_train = [];            
        label_test = [];    
        Y_train_k = [];  
        Y_test_k = [];
        cTest = [];
        cTrain = [];
        c=[];                     

    end


end

% figure
% boxplot(1-enoisy);
figure
boxplot(acc);



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
