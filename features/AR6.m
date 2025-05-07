function [ ar ] = AR6(sample)
% [ ar4_2, ar4_3, ar4_4 ] = AR4(sample)
% Segundo,tereciro,quarto e quinto termo do polinomio de coeficentes autoregressivos de ordem 4.
% sample_fft: FFT das amostras.
% fs: frequência de sample (sample/segundo)

    AR=arburg(sample,6); %AR de ordem 4
    ar6_2=AR(2); %segundo termo da série, pois o primeiro é sempre 1.
    ar6_3=AR(3);
    ar6_4=AR(4);
    ar6_5=AR(5);
    ar6_6=AR(6);
    ar6_7=AR(7);
    
    ar=[ar6_2, ar6_3, ar6_4, ar6_5, ar6_6, ar6_7];

        
end

