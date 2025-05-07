% Daniel Prado de Campos - UTFPR - 20/02/2019
% Processamento de EMG do bovino para classificação entre
% ruminação/ócio/alimentação
% Parte 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Etapas
% 0) Carregamento dos dados pré-definições
% 1) Normalização dos dados: mediana dos picos por animal
% 2) Definição das labels: rotulagem usando vídeo como referência
% 3) Segmentação DTOS e extração de características
% 4) Classificação

% *) Funções

%% 0) Carregamento dos dados e pré-definições
% 

close all
clear 
clc

load DATA

fsa = 2000; %frequência de amostragem
fLa = 10;
fHa = 500;
Ti = 4000;

% abrir dados e definir vetores de tempo
[EMGf1,ta1] = getEMG(EMG1,fsa,fLa,fHa,Ti);
[EMGf2,ta2] = getEMG(EMG2,fsa,fLa,fHa,Ti);
[EMGf3,ta3] = getEMG(EMG3,fsa,fLa,fHa,Ti);
%excluir outliers
EMGf3(6650*fsa:6650.1*fsa)=0;
EMGf3(165.8*fsa:165.9*fsa)=0;
EMGf3(end-20*fsa:end)=[];
ta3(end-20*fsa:end)=[];

%% 1) Normalização
% 

MPD = 0.5; % Distância mínima entre picos: 500 ms 
% não identificar picos menores que 10% do valor máximo
MPH1 = 0.10*max(abs(EMGf1));
MPH2 = 0.10*max(abs(EMGf2));
MPH3 = 0.10*max(abs(EMGf3));

[EMGp1,locs1] = findpeaks(abs(EMGf1),ta1,...
     'MinPeakDistance',MPD,'MinPeakHeight',MPH1);
[EMGp2,locs2] = findpeaks(abs(EMGf2),ta2,...
    'MinPeakDistance',MPD,'MinPeakHeight',MPH2);
 [EMGp3,locs3] = findpeaks(abs(EMGf3),ta3,...
     'MinPeakDistance',MPD,'MinPeakHeight',MPH3);

% normalização 
EMGn1 = EMGf1/mean(EMGp1);
EMGn2 = EMGf2/mean(EMGp2);
EMGn3 = EMGf3/mean(EMGp3);

%% 2) Rotulagem
%

% faltou rotular o dado 2 e 3
% Vou usar os dados do animal 2 apenas, pois foi onde ficou perfeito o
% sinal
% 

%Dimensões das janelas de suavização
Ws = 0.2; 
Wb = 2;
Wf = 5;

EMGn = [EMGn1 EMGn2 EMGn3];
ta = [ta1 ta2+ta1(end)+1 ta3+ta2(end)+ta1(end)+1];

%Suavizar (2 vezes)
EMGsmooth =movmean(abs(EMGn),[Ws*fsa],'includenan');
EMGactiv = movmean(abs(EMGsmooth),[Wb*fsa Wf*fsa],'includenan');
%normalização
EMGactiv = (EMGactiv-min(EMGactiv))/(max(EMGactiv)-min(EMGactiv)) ;

% plot(ta,EMGn)
% hold
% plot(ta,EMGsmooth)
% plot(ta,EMGactiv)
% %Se for maior que 20%, é sinal
% plot(ta,EMGactiv>0.2)

%Os labels não são 100%
%labels:
N = length(EMGn); 
T = zeros(N,3); %Target
Tchew = EMGactiv>0.2; %Se >20%, é sinal

N1 = length(EMGn1);
N2 = length(EMGn2);
N3 = length(EMGn3);

t1eat = [1:1100 1975:2600 3050:3350 4200:4415]*fsa;
t1rum = [5600:5914];

t2eat = [1:2200]*fsa+N1;
t2rum = [5000:6100]*fsa+N1;

t3eat = [32:450 550:660 1250:1520 1900:2210 ...
    2650:3000 3650:3820 6040:6250]*fsa + N1 + N2;
t3rum = [5130:6015]*fsa + N1 + N2;

T([t1eat t2eat t3eat],1)=Tchew([t1eat t2eat t3eat])'; %Período de alimentação (visual)
T([t1rum t2rum t3rum],2)=Tchew([t1rum t2rum t3rum])'; %Período de ruminação
T(:,3)=~Tchew; %O que nao for um dos dois, é ócio




%% 3) Segmentação DTOS e extração de características

%Definição dos tamanhos de janela e quantidade de sobreposição
Wi = [0.6 0.8 1 2 3]*fsa; %400 ms ~3 s
ki = [1 2 3]; %100%/50%/25%
J = length(Wi);
K = length(ki);
Wcrit = 0.1*fsa;
%fHtke = 50;
win = 0.05*fsa;
BL = [1.7 1.8]*1e7;

[~,ex]=energyop(EMGn',false);
x = smooth(abs(ex),win,'lowess');

for j = 1:J
    for k=1:K
tic
%Calcula dimensão da janela
W = Wi(j);


% Usar pra obter os segmentos para extração posterior 
%[segT{j,k},wT{j,k},FT{j,k}] = DTSeg(EMGn', W, Wcrit, k,T,fHtke,fsa,BL);
[wT{j,k},FT{j,k}] = DTSeg(EMGn', W, Wcrit,ki(k),T,x,BL);

 toc



    end
end





%% 4) Classificação
%
I = eye(3);

for j = 1:J
    for k=1:K

%Testar com conjunto de Hudgins
YT= wT{j,k};
XT = FT{j,k}(1:4,:)'; 
Tmat = I(YT,:);
nI = min(sum(Tmat));
YTb = [];
XTb = [];
for i = 1:3
    
    Ti = find(Tmat(:,i));
    R = randperm(length(Ti),nI);
    YTb = cat(1,YTb,YT(Ti(R))');
    XTb = cat(1,XTb,XT(Ti(R),:));
    
    Ti = [];
    R = [];
end



%Testar com LDA
mdl_LDA_T = fitcdiscr(XTb,YTb);
e(:,j,k) = kfoldLoss(crossval(mdl_LDA_T),'Mode','individual');

    end
end




%% Plotar resultados
%
a = 0.25; %distância entre box
bW = 0.15; %largura do box
position = [];
tickspos=[];
color = ['y', 'm', 'c']; %cor do box
C=[];
for i = 1:J
   position = cat(2,position,[i-a i i+a]); 
   tickspos = cat(2,tickspos,mean(position(1+(i-1)*K:3+(i-1)*K))); 
   C = cat(2,C,color);
end

figure
boxplot([(1-e(:,:,1)) (1-e(:,:,2)) (1-e(:,:,3))],...
    [ K*(1:J)-2 K*(1:J)-1 K*(1:J)],'widths',bW,'positions',position,'Colors','kkk')

set(gca,'xtick',tickspos)
set(gca,'xticklabel',{'600 ms','800 ms','1 s','2 s','3 s','5s','10s'})

h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),C(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:3), '5 \sigma','10 \sigma','15 \sigma');

set(gca,'xtick',tickspos)
set(gca,'xticklabel',{'600 ms','800 ms','1 s','2 s','3 s','5s','10s'})



%% Funções
% 

function [EMGf,t]=getEMG(EMG,fs,fL,fH,Ti)
EMGf = filtroEMG(EMG,fs,fH,fL);
EMGf(1:Ti)=0;
t = linspace(0,length(EMGf)/fs,length(EMGf));
end

function [segM,w,F] = BlindSeg(data,W,s,T)
% Blind segmentation (fixed length)
% Segment each step.
%
% Inputs:
% - data: data vector
% - W: Window size
% - s: Step size 
% Outputs:
% - segM: segment matrix

N = length(data);
segM = zeros(round(N/s)+2,W);

k=1;
w = [];
for i = 1:s:N-W
segM(k,:) = data(i:i+W-1);

Tw = T(i:i+W,:);
[~,wi]=max(sum(Tw));
w=cat(2,w,wi);

th_ssc=0.0003;
th_zc=0.0001;
%fs = 2000;

F(1,k)=MAV(segM);
F(2,k)=WL(segM);
F(3,k)=ZC(th_zc,segM);
F(4,k)=SSC(th_ssc,segM);

k=k+1;
end

end


function [wT,F] = DTSeg(data, W, Wcrit, k,T,x,BL)
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
th = k*std(x(BL(1):BL(2)));

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
   if (C>Wcrit) && (Ci+W-1)<N
       %We got a new segment!
       w=Ci:(Ci+W-1);             %segment boundary
       
       Tw = T(w,:);
       [~,wi]=max(sum(Tw));
       wT=cat(2,wT,wi);
       
       %where(w) = 1;        %If there is a segment, there is a 1.
       %segM(j,:) = data(w);   %Add the segment to the cell
       segM = data(w);   %Add the segment to the cell
       
       
th_ssc=0.0003;
th_zc=0.0001;
%fs = 2000;

%Hudgins: MAV, WL, ZC, SSC
%TD4: LS, MFL, MSR, WAMP
%TD9: LS, MFL, MSR, WAMP, ZC, RMS, IAV, DASDV, VAR


F(1,j)=MAV(segM);
F(2,j)=WL(segM);
F(3,j)=ZC(th_zc,segM);
F(4,j)=SSC(th_ssc,segM);


       j=j+1;
       
       
    
   end
   
end

end

