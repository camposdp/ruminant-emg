function [ wamp ] = WAMP(wamp_th,sample)
% WAMP(wamp_th,sample)
% Esse atributo é utilizado para medir a amplitude Willison.
% Está relacionado com a mudança brusca de valores.
% wamp_th: É o critério de amplitude.
% sample: É o vetor da amostra.

N=length(sample)-1;
wamp=0;

for k=2:N
    
   if abs(sample(k)-sample(k-1)) > wamp_th
    
         wamp = wamp + 1;
    end
        
        
end

