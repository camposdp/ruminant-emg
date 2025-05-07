% Daniel Prado de Campos - UTFPR - 22/04/2019
% Preparar validação.teste e k-fold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



close all
clear 
clc

load EMG_raw_labels
load EMGre_exp

fsa = 2000; %frequência de amostragem

BL = [4.45 4.7]*1e5;


W = [0 0.6]*fsa;  

Wcrit = [0.10 1]*fsa;

Nb = W(1);
Nf = W(2);

%% 3) Segmentação DTOS e extração de características


Wsmooth = 0.15*fsa;
Kstd = 4.0;


x = smoothdata(abs(EMGn'),'movmean',Wsmooth);
th = mean(EMGt_re(BL(1):BL(2)))+Kstd*std(EMGt_re(BL(1):BL(2)));


% Usar pra obter os segmentos para extração posterior 
[xi,xf,wT] = DTSeg(EMGn',Wcrit,x,th,T);

xi = xi(1:end-1);
xf = xf(1:end-1);
wT = wT(1:end-1);



Clab{1} = find(wT==1);
Clab{2} = find(wT==2);
Clab{3} = find(wT==3);
K = 0.8;

    wTv = [];
    xiv = [];
    wTt = [];
    xit = [];

for i = 1:2
    %N = length(Clab{i});
    N = 1500;
    Nval = randperm(N,round(N*K));
    Ntest = 1:N; 
    Ntest(Nval) = [];
    
    wTv = [wTv wT(Clab{i}(Nval))];
    xiv = [xiv xi(Clab{i}(Nval))];
    wTt = [wTt wT(Clab{i}(Ntest))];
    xit = [xit xi(Clab{i}(Ntest))];
end


label_test = wTt;
label_train = wTv;


N_t = length(label_test);
N_v = length(label_train);

%% 1) Preparar Matriz de EMG

w0 = W(1)+W(2);

d=1;
w = w0/d;
Y_train = zeros(w,N_v);
Y_test = zeros(w,N_t);

for i = 1:N_v
Y_train(:,i) = decimate(EMGn(xiv(i)-Nb:xiv(i)+Nf-1),d); 
end

for i = 1:N_t
Y_test(:,i) = decimate(EMGn(xit(i)-Nb:xit(i)+Nf-1),d); 
end


Ck = cvpartition(N_v,'KFold',10);


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
