function [ conta ] = ZC(amp_th,sample)
% ZC(amp_th,sample)
% Esse atributo � utilizado para detectar a passagem por zero:
% amp_th: � o crit�rio para discrminar o ru�do.
% sample: � o vetor da amostra. 
N=length(sample)-1;
conta=0;

for k=2:N
    
     Di=sample(k-1)-sample(k);
     Fi=sample(k-1)*sample(k);
     S=-sign(Fi);
     
     if (Di > amp_th) && (S>0)
         
         conta = conta + 1;
     
     end
        
end
        
        
end


