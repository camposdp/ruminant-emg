close all
clear 
clc

load EMG_raw_labels
load EMGre_exp
load Segments_labels_Cpartition

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

for i = 1:length(wT)

Y = EMGn(xi(i)-Nb:xi(i)+Nf-1); 
for j = 1:2400
e1(i,j) = abs(sum(Y_train(:,j) - Y')); 
end

for j = 1:600
e2(i,j) = abs(sum(Y_test(:,j) - Y')); 
end

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