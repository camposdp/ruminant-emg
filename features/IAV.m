function [ iav ] = IAV(sample)
% MAV(sample)
% Esse atributo � utilizado para medir o valor absoluto m�dio.
% sample: � o vetor da amostra. 

N=length(sample);
iav=0;
for k=1:N
    
   iav = iav + abs(sample(k)); 
        
end
        
        
end

