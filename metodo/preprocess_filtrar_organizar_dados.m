

%% Abrir dados, filtrar e normalizar

clear
close all
clc

load('bovino_preproc_sem57_sem74.mat')
load DATA_raw

%Coluna 6
%Linha 17: Silagem (2428) 
%Linha 38: Silagem (2390) 
%Linha 56: Silagem (2000) *RUIM
%Linha 39: Ruminação
%Linha 40: Ruminação
%Linha 57: Ruminação


fsa = 1000; %frequência de amostragem
fLa = 10;
fHa = 499;
Ta = 1000;

[EMG_A1,ta1] = getEMG(D{17,6},fsa,fLa,fHa,Ta,[]); %SIL
[EMG_A2,ta2] = getEMG(D{38,6},fsa,fLa,fHa,Ta,[]); %SIL
[EMG_A3,ta3] = getEMG(D{39,6},fsa,fLa,fHa,Ta,[]); %RUM
[EMG_A4,ta4] = getEMG(D{40,6},fsa,fLa,fHa,Ta,[]); %RUM
[EMG_A5,ta5] = getEMG(D{57,6},fsa,fLa,fHa,Ta,[]); %RUM


fsb = 2000; %frequência de amostragem
fLb = 10;
fHb = 500;
Tb = 4000;



%excluir outliers

OUT3 = [6650*fsa:6650.1*fsa 165.8*fsa:165.9*fsa];

% abrir dados e definir vetores de tempo
%[EMG_N1,tn1] = getEMG(EMG1,fsb,fLb,fHb,Tb,[]); %ANIMAL PODRE
[EMG_N2,tn2] = getEMG(EMG2',fsb,fLb,fHb,Tb,[]);
[EMG_N3,tn3] = getEMG(EMG3',fsb,fLb,fHb,Tb,OUT3);

EMG_N3(end-20*fsa:end)=[];
tn3(end-20*fsa:end)=[];

EMG_N2(1.11165e7:1.114e7)=0;

EMG_N3(3.313e5:3.319e5)=0;
% 
% ANIMAL 1
% 1.345e6 %:silagem
% 1.1218e7:1.1708e7 %:ruminação

% ANIMAL 2
% 2.85e6 %silagem
% 1.075e7:1.209e7 %ruminação

% 1.084

% ANIMAL 3
% 6.6e4:7.3e5 %silagem
% 2.559e6:2.797e6 %silagem
% 1.024e7:1.203e7 %ruminação


 
ta6 = tn2(1:2.85e6);               %SIL
ta7 = tn2(1.075e7:1.209e7);      %RUM

ta8 = tn3(2.559e6:2.797e6);        %SIL
ta9 = tn3(1.024e7:1.203e7);       %RUM

ta10 = tn3(6.6e4:7.3e5);            %SIL

EMG_A6 = EMG_N2(1:2.85e6);               %SIL
EMG_A7 = EMG_N2(1.075e7:1.209e7);      %RUM
EMG_A8 = EMG_N3(2.559e6:2.797e6);        %SIL
EMG_A9 = EMG_N3(1.024e7:1.203e7);       %RUM

EMG_A10 = EMG_N3(6.6e4:7.3e5);       %SIL
EMG_A11 = EMG_N2(3.92e6:4.28e6);       %SIL
EMG_A12 = EMG_N3(5.357e6:5.721e6);       %SIL
EMG_A13 = EMG_N3(3.869e6:4.142e6);       %SIL
EMG_A14 = EMG_N3(7.335e6:7.463e6);       %SIL
%% Funções

function [EMGn,t]=getEMG(EMG,fs,fL,fH,Ti,OUT)
EMGf = filtroEMG(EMG,fs,fH,fL);
EMGf(1:Ti)=0;

EMGf(OUT)=0;

MPD = 0.5; 
MPH = 0.1*max(abs(EMGf));

t = linspace(0,length(EMGf)/fs,length(EMGf));
[EMGp,~] = findpeaks(abs(EMGf),t,...
     'MinPeakDistance',MPD,'MinPeakHeight',MPH);
Y = quantile(EMGp,0.95);
EMGn = EMGf/Y;

end






