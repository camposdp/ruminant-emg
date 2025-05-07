%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daniel P. Campos - UTFPR/CPGEI - 30/05/2019
% Organizar resultados e plotar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Abrir Dados
clear
close all
clc


load Resultados_Handcraft_DadosNovos
load('Resultados_FDDL_DadosNovos_varDic_W_100_400_limiar_novo_KFOLD1000','acc')
acc1 = acc;
load('Resultados_FDDL_DadosNovos_D_500_W_100_400_limiar_novo_KFOLD1000','acc')
acc2 = acc;
%Treino com ruído
%load('Resultados_FDDL_DadosNovos_Noisy','acc')
%Treino sem ruído
load('Resultados_FDDL_DadosNovos_Noisy_Treino_Limpo','acc')
accN = acc;
%Treino com ruído
%load('Resultados_Handcraft_DadosNovos_Noisy')
%Treino sem ruído
load('Resultados_Handcraft_DadosNovos_Noisy_Treino_Limpo')

figure
method = {'MS1','MS2','TD4','TD9','FDDL'};
ACC = 100*[1-ef acc1(5,:)'];

% 
% 
% b=bar(mean(ACC),0.9)
% set(gca, 'XTickLabel',method, 'XTick',1:numel(method))
% hold
% errorbar((1:5),mean(ACC),std(ACC),'k.','LineWidth',1.5)
% ylim([70 95])
% 
% b(1).FaceColor = [.8 .8 .8];

leg = {'MS1','MS2','TD4','TD9','FDDL'};
colors = {'g', 'm', 'c', 'b', 'r'};

h=boxplot(ACC,leg,'color','k','whisker',10)
ylim([65 95])
ax = gca;
ax.YGrid = 'on';

ylabel('Accuracy (%)')


h = findobj(gca,'Tag','Box');

f = flip(1:length(h));
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),...
       colors{f(j)},'FaceAlpha',.5);
end


figure
ACdic = [acc1;acc2];
Dic = [50 100 200 300 400 500];
errorbar(Dic,mean(100*ACdic'),std(100*ACdic'),...
    '-ok','LineWidth',1,'color','r',...
    'MarkerEdgeColor','k','MarkerFaceColor','r')
%hold
%plot(Dic,max(100*AC'),'r--x')
%plot(Dic,min(100*AC'),'b--x')
xticks(Dic);
ylim([78 92])
xlim([0 550])
ylabel('Accuracy (%)')
xlabel('Dictionary Size')
grid on

comb = combnk((1:5),2);
pcrit = 0.05;
for i = 1:length(comb)
    
    p(i) = ranksum(ACC(:,comb(i,1)),ACC(:,comb(i,2)),pcrit);
    d(i) = mean(ACC(:,comb(i,1:2)))/std(ACC(:,comb(i,1:2)));
end

tab = [comb p'];

figure
SNR = (0:20);
acMean = 100*[mean(1-ef1');mean(1-ef2');mean(1-ef3');mean(1-ef4')...
    ;mean(accN')]';



symbols = {'s', 'd', '^', 'v', 'o'};

for i = 1:5
this_spec = ['-', colors{i}, symbols{i}];
plot(SNR, acMean(:,i), this_spec,...
    'MarkerFaceColor',colors{i},...
    'LineWidth',1,...
    'MarkerEdgeColor','k')    
hold on
end

legend(leg,'Location','SouthEast')
xticks(SNR)
xlabel('Signal-to-Noise Rate (SNR)')
ylabel('Accuracy (%)')
ax = gca;
ax.YGrid = 'on';





