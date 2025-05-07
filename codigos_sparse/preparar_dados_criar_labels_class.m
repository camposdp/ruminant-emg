% Daniel 
% Preparar labels
%%
close all
clear 
clc

load DATA_raw
load EMGre_exp

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

t1eat = [1*fsa:1100*fsa 1975*fsa:2600*fsa ...
    3050*fsa:3350*fsa 4200*fsa:4415*fsa];
t1rum = [5600*fsa:5914*fsa];

t2eat = [1*fsa:2200*fsa]+N1;
t2rum = [5000*fsa:6100*fsa]+N1;

t3eat = [32*fsa:450*fsa 550*fsa:660*fsa ...
    1250*fsa:1520*fsa 1900*fsa:2210*fsa 2650*fsa:3000*fsa...
    3650*fsa:3820*fsa 6040*fsa:6250*fsa] + N1 + N2;
t3rum = [5130*fsa:6015*fsa] + N1 + N2;

%T([t1eat t2eat t3eat],1)=Tchew([t1eat t2eat t3eat])'; %Período de alimentação (visual)
%T([t1rum t2rum t3rum],2)=Tchew([t1rum t2rum t3rum])'; %Período de ruminação
%T(:,3)=~Tchew; %O que nao for um dos dois, é ócio

T([t1eat t2eat t3eat],1) = 1;
T([t1rum t2rum t3rum],2) = 1;
Tact = [t1eat t2eat t3eat t1rum t2rum t3rum];
Tidle = 1:(N1+N2+N3);
Tidle(Tact)=[];
T(Tidle,3)=1;






function [EMGf,t]=getEMG(EMG,fs,fL,fH,Ti)
EMGf = filtroEMG(EMG,fs,fH,fL);
EMGf(1:Ti)=0;
t = linspace(0,length(EMGf)/fs,length(EMGf));
end

