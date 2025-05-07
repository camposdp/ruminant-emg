function [ ar ] = AR9(sample)
% [ ar4_2, ar4_3, ar4_4 ] = AR4(sample)
% Segundo,tereciro,quarto e quinto termo do polinomio de coeficentes autoregressivos de ordem 4.
% sample_fft: FFT das amostras.
% fs: frequência de sample (sample/segundo)

    AR=arburg(sample,9); %AR de ordem 4
    ar9_2=AR(2); %segundo termo da série, pois o primeiro é sempre 1.
    ar9_3=AR(3);
    ar9_4=AR(4);
    ar9_5=AR(5);
    ar9_6=AR(6);
    ar9_7=AR(7);
    ar9_8=AR(8);
    ar9_9=AR(9);
    
    ar=[ar9_2, ar9_3, ar9_4, ar9_5, ar9_6, ar9_7 ar9_8 ar9_9];

        
end