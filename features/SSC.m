function [ conta ] = SSC(slope_th,sample)
% SSC(slop_th,sample)
% Esse atributo � utilizado para detectar quantos vezes existe a varia��o
% da inclina��o do sinal:
% slope_th: � o crit�rio da discrminar pequenas varia��es.
% sample: � o vetor da amostra. 
%N=length(sample)-1;
%conta=0;

f=(sample(2:end-1)-sample(1:end-2)).*(sample(2:end-1)-sample(3:end));

conta=sum(f>slope_th);
% 
% for k=2:N-1
%     
%     f=(sample(k)-sample(k-1))*(sample(k)-sample(k+1));
%         
%         if f > slope_th;
%         
%             conta = conta + 1;
%             
%         end
%         
% end
%         
        
end

