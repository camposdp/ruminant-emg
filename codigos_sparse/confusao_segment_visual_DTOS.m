
% Daniel Prado de Campos - UTFPR - 22/04/2019
% Avaliação da segmentação DTOS - Confusão
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%


load EMG_filt_e_norm
load EMG_visual
%load EMG_segment_labels
load cont



Xi = xi;
Xf = xf;
X = (Xf - Xi)/fsa*1000;
EMGt_re = EMGt(1:xf(end));




fsa = 2000;
Wcrit = [0.015*fsa 0.6*fsa];
BL = [4.45 4.7]*1e5;


%[Yi,Yf] = DTSeg(EMGt_re, Wcrit,x,th);
%Y = (Yf - Yi)/fsa*1000;
%[P]=getGausMix(Y,St);


Xbin = zeros(length(EMGt_re),1);

for i = 1:length(X)
    Xbin(round(Xi(i)):round(Xf(i)))=1;
end
Xbin = Xbin(1:length(EMGt_re));

Kvar = [0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.5 3 3.5 4 5 6 7 8 9 10];
Wvar = [0.01 0.02 0.03 0.05 0.1 0.15];
for j = 1:length(Wvar) 
for l = 1:length(Kvar)
    tic
    W = Wvar(j)*fsa;
    x = smoothdata(abs(EMGt_re),'movmean',W);
    k=Kvar(l);
    th = mean(x(BL(1):BL(2)))+k*std(x(BL(1):BL(2)));
    [Yi,Yf] = DTSeg(EMGt_re, Wcrit,x,th);
    Y = (Yf - Yi)/fsa*1000;
    Ybin = zeros(length(EMGt_re),1);

    for i = 1:length(Y)
        Ybin(Yi(i):Yf(i))=1;
    end

    [c(l,j),~,~,~] = confusion(Xbin',Ybin');
        toc
end
end

figure
plot(Kvar,100*c,'x-');
xlabel('k (threshold = \mu + k*\sigma)')
ylabel('Error (%)')
legendCell = cellstr(num2str(Wvar', 'W=%1.2f'));
legend(legendCell,'Location','SouthEast')
xticks(Kvar)

[cmin,argmin] = min(c);
axes('Position',[.6 .68 .3 .2])
box on
plot(Wvar,cmin*100,'-or','LineWidth',1)
title('minimum Error')
xlabel('Smooth Window Size')
ylabel('Error (%)')
xlim([0 0.17])
ylim([5 7])
xticks(Wvar)
xtickangle(45)
tickCell = cellstr(num2str(Kvar(argmin)', 'k=%1.1f'));
xt = Wvar;
xt(3:end) = xt(3:end)+0.003;
yt = cmin*100;
yt(3:end) = yt(3:end)-0.2;
yt(1:2) = yt(1:2)+0.2;
text(xt',yt',tickCell);



function [xi,xf] = DTSeg(data,Wcrit,x,th)
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
%wT = [];

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
       %w=Ci:(Ci+W-1);             %segment boundary
       
       %Tw = T(w,:);
       %[~,wi]=max(sum(Tw));
       %wT=cat(2,wT,wi);
       
       %where(w) = 1;        %If there is a segment, there is a 1.
       %segM(j,:) = data(w);   %Add the segment to the cell
       %segM = data(w);   %Add the segment to the cell
       
       
       
       j=j+1;
       
       
    
   end
   
end

end



