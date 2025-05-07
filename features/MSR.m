function [msr] = MSR(x)
%Mean value of the Square Root

N = length(x);
soma = 0;

for i = 1:N
    soma = soma + sqrt(abs(x(i)));
    
end
msr = soma/N;

end

