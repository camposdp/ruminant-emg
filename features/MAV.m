function [ mav ] = MAV(sample)
% MAV(sample)
% Esse atributo é utilizado para medir o valor absoluto médio.
% sample: É o vetor da amostra. 

N=length(sample);
mav=0;
for k=1:N
    
   mav = mav + abs(sample(k))/N; 
        
end
        
        
end

