function [ comprimento ] = WL(sample)
% WL(sample)
% author: Daniel P Campos - UTFPR
% Esse atributo � utilizado para medir o comprimento da curva:
% sample: � o vetor da amostra. 

N=length(sample)-1;
comprimento=0;


for k=2:N
    
  comprimento = comprimento + abs(sample(k)-sample(k-1)) ;
        
end
        
        
end

