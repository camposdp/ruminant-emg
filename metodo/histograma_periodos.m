%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daniel P. Campos - UTFPR/CPGEI - 24/05/2019
% Segmentar e tirar estatísticas de períodos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Abrir Dados
clear
close all
clc


load EMG_base_exp
load EMGre_exp
load KFOLD_1000 

%load KFOLD_500

%% Definições

fsa = 2000;

BL = [4.45 4.7]*1e5;

Wcrit = [0.1*fsa 1*fsa];
Wsmooth = 0.03*fsa;
k=2.5;


EMG_R = [EMG_A7;EMG_A9];
EMG_S = [EMG_A6;EMG_A8;EMG_A10;EMG_A11;EMG_A12;EMG_A13;EMG_A14];

%% Segmentação

x_r = smoothdata(abs(EMG_R),'movmean',Wsmooth);
x_s = smoothdata(abs(EMG_S),'movmean',Wsmooth);

th = mean(EMGt_re(BL(1):BL(2)))+k*std(EMGt_re(BL(1):BL(2)));

[xiR,xfR] = DTOSSeg(EMG_R,Wcrit,x_r,th);
[xiS,xfS] = DTOSSeg(EMG_S,Wcrit,x_s,th);

% Extrair estatísticas

TR = (1e3/fsa)*(xfR - xiR);
TS = (1e3/fsa)*(xfS - xiS);

TiR = (1e3/fsa)*(xiR(2:end)-xfR(1:end-1));
TiS = (1e3/fsa)*(xiS(2:end)-xfS(1:end-1));

BW = 12;
%BWi = 30;
BWi = 500

figure
hold all
histogram(TR,'BinWidth',BW,'Normalization','pdf')
histogram(TS,'BinWidth',BW,'Normalization','pdf')
xlabel('Período (ms)')
ylabel('Densidade de Probabilidade')
legend('Ruminação','Ingestão')
xlim([50 350])

figure
hold all
histogram(TiR,'BinWidth',BWi,'Normalization','pdf')
histogram(TiS,'BinWidth',BWi,'Normalization','pdf')
xlim([2000 10000])
xlabel('Intervalo (ms)')
ylabel('Densidade de Probabilidade')
legend('Ruminação','Ingestão')



function [xi,xf] = DTOSSeg(data,Wcrit,x,th)
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
%       w=Ci:Cf;             %segment boundary
       %w_fix=Ci:(Ci+W-1)
       
%        Tw = T(w,:);
%        [~,wi]=max(sum(Tw));
%        wT=cat(2,wT,wi);
       
       %where(w) = 1;        %If there is a segment, there is a 1.
       %segM(j,:) = data(w);   %Add the segment to the cell
       %segM = data(w);   %Add the segment to the cell
       
       
       
       j=j+1;
       
       
    
   end
   
end

end










