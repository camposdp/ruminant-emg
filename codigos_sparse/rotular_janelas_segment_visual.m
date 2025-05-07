% Daniel Prado de Campos - UTFPR - 22/04/2019
% Rotulagem manual das janelas de segmentação
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Etapas
% 0) Carregamento dos dados pré-definições
% 1) Contagem

% *) Funções
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0) Carregamento dos dados pré-definições


close all
clear 
clc

load EMG_filt_e_norm

fsa = 2000; %frequência de amostragem

T = [3.756e5:5.151e5 ...
    1.122e7:1.131e7 ...
    1.29e7:1.38e7 ...
     2.277e7:2.295e7 ...
     3.09e7:3.092e7 ...
     3.644e7:3.663e7...
     3.734e7:3.743e7...
     1.2e7:1.204e7...
     2.3825e7:2.393e7...
     3.05e4:7.7e4];

 
 

 %% 1) Contagem
 
figure('units','normalized','outerposition',[0 0 1 1])

 EMGt = EMGn(T);

W = 10e3;

button = 1;
i=1;
k=1;

while k+W < length(EMGt) % && (sum(button) <=1
        
        plot(abs(EMGt))
        ylim([0 2])
        xlim([k k+W])
        
        [x,~,button] = ginput(1) 
        if button == 1 
        xi(i)=x;
        [xf(i),~,button] = ginput(1)
        k = xf(i);
        i = i +1;
        else
            k=k+W;
        
        end      

end


xf(end) = [];
xi(end) = [];

EMGt_re = EMGt(1:xf(end));

figure
hold all
plot(EMGt_re)
stem(xf,ones(length(xf),1),'bo')
stem(xi,ones(length(xi),1),'rx')
legend('Inicio da Janela','Fim da Janela')












