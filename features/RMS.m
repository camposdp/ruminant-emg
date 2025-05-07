function [ rms ] = RMS(sample)
% RMS(sample)
% Esse atributo é utilizado para medir o valor RMS do sinal.
% sample: É o vetor da amostra.

N=length(sample)-1;
rms=0;
ssum=0;

for k=2:N
   
    ssum = ssum + sample(k)^2;
    
        if k==N
          
        rms=sqrt(ssum/N);    
            
        end
        
        
end

