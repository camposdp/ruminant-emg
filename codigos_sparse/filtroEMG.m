function [ CH_S ] = filtroEMG(sample,fs,fL,fH)
%filtrop Filtro padrão para EMG
%author: Daniel P Campos - UTFPR
%   filtrop(sample,fs,fL,fH)
%   Seleciona banda de frequência e filtra ruídos da rede
%   -> sample: amostra (sinal)
%   -> fs: frequência de amostragem
%   -> fL: frequência inferior da banda
%   -> fH: frequência superior da banda

CH_S=sample;
 
%Definição dos parâmetros das funções dos filtros, ver a documentação da 
%função de cada filtro:
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

