% Daniel Prado de Campos - UTFPR - 28/02/2019
% Pós-Processamento de EMG do bovino para classificação entre
% ruminação/ócio/alimentação
% Parte 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Rótulos
%1: Eating
%2: Rumination
%3: Idle


clear
clc

%Carregar F e wi
load DATA_features_and_labels
%load DATA_FFT


close all


WinS = 5; %2s 
Ovl = 1; %disr.

%Primeiras 6 features
f = F{WinS,Ovl}([1:6],:);
wi=w{WinS,Ovl};

%% PCA
%

FeatName = {'MAV','WL','ZC','SSC','SSI','LD'};
Class = {'Eating','Rumination','Idle'};

[coeff,score,latent,~,explained] = pca(f');
expPCA = explained(1) + explained(2);
Xcentered = score*coeff';
figure
hold all
ci = Class(wi);

hbi=biplot(coeff(:,1:2),'scores',score(:,1:2),...
    'varlabels',FeatName,'ObsLabels',ci);
% Manipulate plot colors
col = [1 0 0; 0 1 0; 0 0 1];
dataOffset = length(hbi)-length(wi);
for i = 1:length(wi)

        
            set(hbi(i+dataOffset), 'Color', col(wi(i),:));

        
            set(hbi(i+dataOffset), 'MarkerSize', 3);
            set(hbi(i+dataOffset), 'Marker', 'o');
            set(hbi(i+dataOffset), 'LineWidth', 0.5);
end

xlabel({'PC1: ',num2str(explained(1),'%2.2f'),'%'});
ylabel({'PC2: ',num2str(explained(2),'%2.2f'),'%'});




%Se quiser plotar 3d:

%Quais features?
% a=6;
% b=3;
% c=4;

% figure
% scatter3(f(a,wi==1),f(b,wi==1),f(c,wi==1),'*')
% hold all
% scatter3(f(a,wi==2),f(b,wi==2),f(c,wi==2),'+')
% scatter3(f(a,wi==3),f(b,wi==3),f(c,wi==3),'.')
% xlabel(FeatName(a))
% ylabel(FeatName(b))
% zlabel(FeatName(c))
% legend(Class)


%% Avaliação dos Feature Sets e Classificadores
%

Y = wi;
X = f; 

%Definição dos sets
FeatSets = {'MAV','WL','ZC','SSC','SSI','LD',...
    'ZC+SSC','ZC+SSI','PC1+PC2','Hudgins','All'};
x{1} = X(1,:);
x{2} = X(2,:);
x{3} = X(3,:);
x{4} = X(4,:);
x{5} = X(5,:);
x{6} = X(6,:);
x{7} = X(3:4,:);
x{8} = X([3 5],:);
x{9} = score(:,1:2)';
x{10} = X(1:4,:);
x{11} = X(1:6,:);


for i = 1:length(x)
    Xi = x{i};

 %Cria modelo do classificador
k=10; %k do kNN
mdl_LDA = fitcdiscr(Xi',Y);
mdl_KNN = fitcknn(Xi',Y,'NumNeighbors',k);
mdl_DT = fitctree(Xi',Y);

%Erro de classificação com k-fold cross validation
e(:,i,1) = kfoldLoss(crossval(mdl_LDA),'Mode','individual');
e(:,i,2) = kfoldLoss(crossval(mdl_KNN),'Mode','individual');
e(:,i,3) = kfoldLoss(crossval(mdl_DT),'Mode','individual');

end


%% Plotar resultados
%

s= size(e);
K = s(3);
J = s(2);

a = 0.2; %distância entre box
bW = .10; %largura do box
position = [];
tickspos=[];
color = ['y', 'm', 'c']; %cor do box
C=[];
for i = 1:J
   position = cat(2,position,[i-a i i+a]); 
   tickspos = cat(2,tickspos,mean(position(1+(i-1)*K:3+(i-1)*K))); 
   C = cat(2,C,color);
end

figure
boxplot([(1-e(:,:,1 )) (1-e(:,:,2)) (1-e(:,:,3))],...
    [ K*(1:J)-2 K*(1:J)-1 K*(1:J)],'widths',bW,'positions',position,'Colors','kkk')

set(gca,'xtick',tickspos)
set(gca,'xticklabel',FeatSets)

h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),C(j),'FaceAlpha',.5);
end

ClassName = {'LDA','KNN','DT'};
c = get(gca, 'Children');
hleg1 = legend(c(1:3), ClassName);


%% Teste de Wilcoxon
% 

Cfeat=combnk([1:5],2);
Cclass=combnk([1:3],2);

for j = 1:3
for i = 1:length(Cfeat)
X = e(:,6+Cfeat(i,1),j);
Y = e(:,6+Cfeat(i,2),j);

Pfeat(i,j) = ranksum(X,Y);
end
end

for j = 1:length(Cclass)
for i = 1:5
X = e(:,6+i,Cclass(j,1));
Y = e(:,6+i,Cclass(j,2));

Pclass(i,j) = ranksum(X,Y);
end
end

% Montar tabela do teste
%

T_class = [];
tblVarCls = ClassName(Cclass); 
for i=1:size(Pclass,2)
varName = [tblVarCls{i,1} '_' tblVarCls{i,2}];
T_temp = table(Pclass(:,i),'VariableNames',{varName});
T_class = [T_class T_temp];
end
T_class.Properties.RowNames = FeatSets(7:11);


T_feat = [];
tblVarFeat = FeatSets(Cfeat+6);
for i=1:size(Pfeat,2)
varName = [ClassName{i}];
T_temp = table(Pfeat(:,i),'VariableNames',{varName});
T_feat = [T_feat T_temp];

end

for i = 1:size(Pfeat,1)
    tRowName{i} =  [tblVarFeat{i,1} '_' tblVarFeat{i,2}];
end
T_feat.Properties.RowNames = tRowName;

T_class % Tabela comparando entre classificadores
T_feat  % Tabela comparando entre feature sets
% Se p < 0.05 rejeita-se a hipótese de que as medianas são iguais 

markSigM(permute(e(:,7:11,:),[3,2,1]),0.05)
squeeze(median(1-e))


%% FFT
%


% figure 
% hold all
% c1 = {'k-','r-','b-'};
% c2 = {'k--','r--','b--'};
% for i = 1:3
% 
% FTmean = mean(WinFFT(wi==i,:));
% FTmax = max(WinFFT(wi==i,:));
% %FTmin = min(WinFFT(wi==i,:));
% plot(fi{WinS,Ovl},FTmean,c1{i});
% %plot(fi{WinS,Ovl},FTmax,c2{i});
% %plot(fi{WinS,Ovl},FTmin,'r--');
% legend(Class(i))
% end
% legend(Class)
