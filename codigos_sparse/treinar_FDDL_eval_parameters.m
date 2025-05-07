close all
clear
clc

load Segments_labels_Cpartition
load Resultados_FDDL_varDic_W_0_400
load EMG_raw_labels
load EMGre_exp
load Segments_labels_index


fsa = 2000; %frequência de amostragem

%% 3) Segmentação DTOS e extração de características

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




%% 1) Preparar Matriz de EMG


%Extrair features
idx = find(ind_train==1);
%Usar só wTv e Y_train (o test será usado a posteriori)
wTx = wT(idx);
xix = xi(idx);
YT = wTx;





W = [0.0 0.4]*fsa;

Y_train = [];
for ind = 1:length(idx)
Y_train(:,ind) = EMGn(xix(ind)-W(1):xix(ind)+W(2)-1);
end

    


Y_n = Y_train;


vg = 0.0035;



for i = 1:10



cTest =  Ck.test(i);
label_test = wTx(cTest);    
Y_test = Y_n(:,cTest);    

% cTrain = Ck.training(i);
% label_train = wTx(cTrain);
% Y_train_k = Y_n(:,cTrain);




optsi = opts{i,10};


Di = D{i,10};
CoefMi= CoefM{i,10};

optsi.gamma = vg;
pred = FDDL_pred(Y_test, Di, CoefMi, optsi);
acc = double(numel(find(pred == label_test)))/...
    numel(label_test)
AC(i)=acc;

end




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