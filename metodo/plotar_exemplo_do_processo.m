close all
clc
clear


load 'Resultados_FDDL_DadosNovos_D_500_W_100_400_limiar_novo_KFOLD1000'
load 'Segmented_Signal'

Y = SS(:,300);
[pred,X,r1,r2,E,yr] = FDDL_pred(Y, D{1}, CoefM{1}, opts{1});

t = linspace(0,500,1000);


figure
plot(t,Y,'Color',[0.7 0.7 0.7],'LineWidth',2)
ylabel('sEMG Amplitude')
xlabel('Time (ms)')


figure
subplot(1,2,1)
X1=X(1:500);
A1 = 1:500;
stem(A1(abs(X1)>0.02),X1(abs(X1)>0.02),...
    'b','filled','MarkerSize',3)
title('Eating')
ylabel('Mixture Coefficients')
xlabel('Atoms')
xlim([1 500])

subplot(1,2,2)

X2=X(501:1000);
A2 =501:1000;

stem(A2(abs(X2)>0.02),X2(abs(X2)>0.02),...
    'r','filled','MarkerSize',3)

title('Rumination')
ylabel('Mixture Coefficients')
xlabel('Atoms')
ylim([-0.4 0.4])
xlim([501 1000])

figure
subplot(1,2,1)
hold all
plot(t,Y,'Color',[0.7 0.7 0.7],'LineWidth',2)
plot(t,yr{1},'b')

legend('Signal Sample','Reconstruction from Eating dictionary')
ylabel('sEMG Amplitude')
xlabel('Time (ms)')

subplot(1,2,2)
hold all
plot(t,Y,'Color',[0.7 0.7 0.7],'LineWidth',2)
plot(t,yr{2},'r')

legend('Signal Sample','Reconstruction from Rumination dictionary')
ylabel('sEMG Amplitude')
xlabel('Time (ms)')



figure
hold
b=bar(1,E(1));
set(b,'FaceColor','b');
b=bar(2,E(2));
set(b,'FaceColor','r');
set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'Eating','Rumination'})

ylabel('error (e_i)')

function [pred,X,r1,r2,E,yr] = FDDL_pred(Y, D, CoefM, opts) % GC
    vgamma = opts.gamma;
    opts.max_iter = 300;
    [X, ~] = lasso_fista(Y, D, zeros(size(D,2), size(Y,2)), vgamma, opts);
    C = size(CoefM,2);
    % w = 0.5;
    E = zeros(C, size(Y,2));
    for c = 1: C 
        Dc = get_block_col(D, c, opts.D_range);
        Xc = get_block_row(X, c, opts.D_range);
        R1 = Y - Dc*Xc;
        E1 = sum(R1.^2);
        R2 = X - repmat(CoefM(:, c), 1, size(Y,2));
        E2 = sum(R2.^2);
        E(c,:) = E1 + opts.weight*E2;
        r1(c,:)=R1;
        r2(c,:)=R2;
        yr{c} = Dc*Xc; 
        
    end 
    [~, pred] = min(E);
end 