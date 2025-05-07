% Daniel Prado de Campos - UTFPR - 03/05/2019
% Processamento de EMG - Parte 2
% Treinamento kSVD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Etapas
% 0) Carregamento dos dados pré-definições
% 1)Pooling
% 2)Classification


%% 0) Carregamento dos dados pré-definições

clear
clc
close all

load Dic_kSVD
load labels_val_e_test





d = [1 2 4 8];
k = [0.5 1.0 1.5];


for i = 1:length(d)
    for j = 1:length(k)
%     f1 = max(A{j,i});
%     f2 = sum(A{j,i});
    %% 1)Pooling
    X=pooling(A{j,i});
    %F = [sum(abs(A{j,i}))' max(abs(A{j,i}))'];
    %X = full(F);
    Y = wTv;
    %% 2)Classification
     mdl_LDA=fitcdiscr(X,Y);
     e(length(k)*(i-1)+j,:) = kfoldLoss(crossval(mdl_LDA),'Mode','individual');
        
    end
end
boxplot(1-e)


 function F=pooling(As)

 A = full(As);
 k=10;
S = size(A);
N = S(1);
P=permn([1 0],k);
P = P(1:end-1,:);
%P = eye(k);
n = N/k;

for i = 1:k
    ind = [1+n*(i-1):n*i];
    f1(i,:) = sum(abs(A(ind,:)));
    f2(i,:) = max(abs(A(ind,:)));
end

 for q = 1:length(P)
     F1(q,:) = sum(f1(find(P(q,:)),:),1);
     F2(q,:) = max(f2(find(P(q,:)),:),[],1);
 end

F = [F1' F2'];

%F = [f1' f2'];

end



