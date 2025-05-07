close all
%clear
%clc

%abrir EMGn em EMG_filtr_e_norm

t(1,1) = 1.196e7;
t(1,2) = t(1,1) + 90*2000;



t(2,1) = 2.277e7;
t(2,2) = t(2,1) + 90*2000;

t(3,1) = t(1,1) + 21*2000;
t(3,2) = t(3,1) + 3*2000;

t(4,1) = t(2,1)+ 3.1*2000;
t(4,2) = t(4,1) + 3*2000;


[N,~]=size(t);

Y = [3 3 2 2];
X = [90 90 3 3];

for i = 1:N
    figure('Position', [10 10 700 200])
    w = t(i,1):t(i,2);
    T = ta(w) - ta(t(i,1));
    plot(T,EMGn(w)) ;
    ylim([-Y(i) Y(i)])
    xlim([0 X(i)])
    xlabel('Time (s)')
    ylabel('sEMG Amplitude')
    
end

