%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daniel Prado de Campos  - 17/05/2019
% UTFPR-CPGEI
% Comparação de complexidade computacional
% FDDL vs Handcrated features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Carregamento de daddos

close all
clc
clear

load Resultados_FDDL_varDic_W_0_400
%load Resultados_FDDL_varDic_W_0_500
%load Resultados_FDDL_varDic_W_200_400

load EMG_raw_labels
load EMGre_exp
load Segments_labels_Cpartition
load Segments_labels_index


fsa = 2000;
th_ssc=0.001;
th_zc=0.11;
th_wamp = 0.06;

lambda1 = 0.01;
lambda2 = 0.01;    
    
SET = [12 11 10 5 3 6 9 8 7];


%Kdic = 400;
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

W = [0.0 0.4]*fsa;


%% 1) Preparar Matriz de EMG

%Extrair features
idx = find(ind_train==1);

%Usar só wTv e Y_train (o test será usado a posteriori)
wTx = wT(idx);
YT = wTx;



for j = 1:length(idx)
 
id = idx(j);    
%Sv = EMGn(xi(id):xf(id));
Sf = EMGn(xi(id)-W(1):xi(id)+W(2)-1);

Ff(1,j)=ZC(th_zc,Sf);
Ff(2,j)=WAMP(th_wamp,Sf);
Ff(3,j)=RMS(Sf);
Ff(4,j)=VAR(Sf);
Ff(5,j)=DASDV(Sf);
Ff(6,j)=IAV(Sf);
Ff(7,j)=MFL(Sf);
Ff(8,j)=MSR(Sf);
Ff(9,j)=LS(Sf,2);


end

Xtrainf = Ff;






%% 4) Classificação
%
I = eye(3);
Xf = Ff';
mdl_LDA = fitcdiscr(Xf,YT);
i = 1;
%for k = 1:3
for j = 1:length(idx)
 
id = idx(j);       
Sf = EMGn(xi(id)-W(1):xi(id)+W(2)-1);
label_test = wT(id);

Ff = [];
tic
Ff(1)=ZC(th_zc,Sf);
Ff(2)=WAMP(th_wamp,Sf);
Ff(3)=RMS(Sf);
Ff(4)=VAR(Sf);
Ff(5)=DASDV(Sf);
Ff(6)=IAV(Sf);
Ff(7)=MFL(Sf);
Ff(8)=MSR(Sf);
Ff(9)=LS(Sf,2);

loss(mdl_LDA,Ff,label_test);
t1(i) = toc;

tic
FDDL_pred(Sf', D{1,1}, CoefM{1,1}, opts{1,1});
t2(i) = toc;

i = i+1;

end
%end








%% 















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