function [ wamp ] = WAMP(wamp_th,sample)
% WAMP(wamp_th,sample)
% Esse atributo � utilizado para medir a amplitude Willison.
% Est� relacionado com a mudan�a brusca de valores.
% wamp_th: � o crit�rio de amplitude.
% sample: � o vetor da amostra.

N=length(sample)-1;
wamp=0;

for k=2:N
    
   if abs(sample(k)-sample(k-1)) > wamp_th
    
         wamp = wamp + 1;
    end
        
        
end

