function [ ceps ] = CC9(sample)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

   AR=arburg(sample,9); %AR de ordem 4

    c(1)=-AR(1);
    sum=0;
    for i=2:9
       
        for l=1:(i-1)
            
            sum=sum+(1-l/i)*AR(i)*c(i-1);
            
        end    
        
        c(i)=-AR(i)-sum;
        sum=0;
        
    end    
   
    ceps=c(2:end);
    
end

