close all
clc
clear


load 'Resultados_FDDL_DadosNovos_D_500_W_100_400_limiar_novo_KFOLD1000'
load 'Segmented_Signal'


for i = 1:10
Yss = SS(:,i);
Ysr = SR(:,i);
[E1ss(:,i),E2ss(:,i),Ess(:,i)] = FDDL_pred(Yss, D{1}, CoefM{1}, opts{1});
[E1sr(:,i),E2sr(:,i),Esr(:,i)] = FDDL_pred(Ysr, D{1}, CoefM{1}, opts{1});
end


function [E1,E2,E] = FDDL_pred(Y, D, CoefM, opts) % GC
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
        E1(c,:) = sum(R1.^2);
        R2 = X - repmat(CoefM(:, c), 1, size(Y,2));
        E2(c,:) = sum(R2.^2);
        E(c,:) = E1(c,:) + opts.weight*E2(c,:);

        
    end 
    [~, pred] = min(E);
end 