function [ var ] = VAR(sample)
% RMS(sample)
% Esse atributo é utilizado para medir o valor RMS do sinal.
% sample: É o vetor da amostra.

N=length(sample);

sum=0;

for k=1:N
   
    sum = sum + sample(k)^2;
   
        
end
          
        var=sum/(N-1);    
           

end