function [ CH_S ] = filtroEMG(sample,fs,fL,fH)
%filtrop Filtro padr�o para EMG
%author: Daniel P Campos - UTFPR
%   filtrop(sample,fs,fL,fH)
%   Seleciona banda de frequ�ncia e filtra ru�dos da rede
%   -> sample: amostra (sinal)
%   -> fs: frequ�ncia de amostragem
%   -> fL: frequ�ncia inferior da banda
%   -> fH: frequ�ncia superior da banda

CH_S=sample;
 
%Defini��o dos par�metros das fun��es dos filtros, ver a documenta��o da 
%fun��o de cada filtro:
%help butter
%help iirnotch

PA=fL*2/fs;  
PB=fH*2/fs;

Wo=50*2/fs;
Q=30;
BW=Wo/Q;

 [Ba,Aa] = butter(4,PA,'high');  %Filtro Butterworth de ordem 3 passa-alta
  %FPA=tf(Ba,Aa,1/1000); 
 CH_S=filter(Ba,Aa,CH_S); %Aplica o filtro
 
 
 [Bb,Ab] = butter(4,PB,'low');  %Filtro Butterworth de ordem 3 passa-baixa
 % FPB=tf(Bb,Ab,1/1000);
 CH_S=filter(Bb,Ab,CH_S); %Aplica o filtro

 [Bn,An] = iirnotch(Wo,BW); % Filtro "Notch" rejeita-faixa em 60 Hz;
% FN=tf(Bn,An,1/1000); 
 CH_S=filter(Bn,An,CH_S); %Aplica o filtro
 
  [Bn,An] = iirnotch(2*Wo,2*BW); % Filtro "Notch" rejeita-faixa em 60 Hz;
% FN=tf(Bn,An,1/1000); 
 CH_S=filter(Bn,An,CH_S); %Aplica o filtro

   [Bn,An] = iirnotch(3*Wo,BW); % Filtro "Notch" rejeita-faixa em 60 Hz;
% FN=tf(Bn,An,1/1000); 
 CH_S=filter(Bn,An,CH_S); %Aplica o filtro

    [Bn,An] = iirnotch(4*Wo,BW); % Filtro "Notch" rejeita-faixa em 60 Hz;
% FN=tf(Bn,An,1/1000); 
 CH_S=filter(Bn,An,CH_S); %Aplica o filtro
 
    [Bn,An] = iirnotch(5*Wo,BW); % Filtro "Notch" rejeita-faixa em 60 Hz;
% FN=tf(Bn,An,1/1000); 
 CH_S=filter(Bn,An,CH_S); %Aplica o filtro
%  
%      [Bn,An] = iirnotch(6*Wo,BW); % Filtro "Notch" rejeita-faixa em 60 Hz;
% % FN=tf(Bn,An,1/1000); 
%  CH_S=filter(Bn,An,CH_S); %Aplica o filtro
%  
%      [Bn,An] = iirnotch(7*Wo,BW); % Filtro "Notch" rejeita-faixa em 60 Hz;
% % FN=tf(Bn,An,1/1000); 
%  CH_S=filter(Bn,An,CH_S); %Aplica o filtro
%  
%      [Bn,An] = iirnotch(8*Wo,BW); % Filtro "Notch" rejeita-faixa em 60 Hz;
% % FN=tf(Bn,An,1/1000); 
%  CH_S=filter(Bn,An,CH_S); %Aplica o filtro
%  
%      [Bn,An] = iirnotch(9*Wo,BW); % Filtro "Notch" rejeita-faixa em 60 Hz;
% % FN=tf(Bn,An,1/1000); 
 CH_S=filter(Bn,An,CH_S); %Aplica o filtro
 
end

