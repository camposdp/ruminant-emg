function [ mav ] = MAV(sample)
% MAV(sample)
% Esse atributo � utilizado para medir o valor absoluto m�dio.
% sample: � o vetor da amostra. 

N=length(sample);
mav=0;
for k=1:N
    
   mav = mav + abs(sample(k))/N; 
        
end
        
        
end

