function [ LD ] = LD(sample)
% LD(sample)

N=length(sample);
x=0;

for k=1:N
    
 x = x + log(abs(sample(k)));
        
end

LD = exp(x/N);

end
