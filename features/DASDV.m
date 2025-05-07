function [ dasdv ] = DASDV(sample)
% DASDV


N=length(sample)-1;
somaD=0;
for k=1:N-1
    
  somaD = somaD + (sample(k+1)-sample(k))^2 ;
        
end
        dasdv = sqrt(somaD/(N-1));
        
end

