function [ iav ] = IAV(sample)
% MAV(sample)
% Esse atributo é utilizado para medir o valor absoluto médio.
% sample: É o vetor da amostra. 

N=length(sample);
iav=0;
for k=1:N
    
   iav = iav + abs(sample(k)); 
        
end
        
        
end

