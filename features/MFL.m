function [mfl] = MFL(x)
%maximum fractal length 

N = length(x);
soma = 0;

for i = 1:N-1
    soma = soma + (x(i+1)-x(i))^2;
end

mfl = log10(sqrt(soma));


end

