% Daniel Prado de Campos 03/05
% Preparar validação - 80-20
%%%


close all
clear
clc

load EMG_labels_and_windows


C{1} = find(wT==1);
C{2} = find(wT==2);
C{3} = find(wT==3);
K = 0.8;

    wTv = [];
    xiv = [];
    wTt = [];
    xit = [];

for i = 1:3
    N = length(C{i});
    Nval = randperm(N,round(N*K));
    Ntest = 1:N; 
    Ntest(Nval) = [];
    wTv = [wTv wT(C{i}(Nval))];
    xiv = [xiv xi(C{i}(Nval))];
    wTt = [wTt wT(C{i}(Ntest))];
    xit = [xit xi(C{i}(Ntest))];
end