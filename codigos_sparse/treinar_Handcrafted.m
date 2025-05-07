% Daniel Prado de Campos - UTFPR - 20/02/2019
% Processamento de EMG do bovino para classificação entre
% ruminação/ócio/alimentação
% Parte 1 - Hand-crafted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Etapas
% 0) Carregamento dos dados pré-definições
% 1) Normalização dos dados: mediana dos picos por animal
% 2) Definição das labels: rotulagem usando vídeo como referência
% 3) Segmentação DTOS e extração de características


% *) Funções

%% 0) Carregamento dos dados e pré-definições
% 

close all
clear 
clc

load EMG_raw_labels
load EMGre_exp
load Segments_labels_Cpartition
load Segments_labels_index


fsa = 2000; %frequência de amostragem

%% 3) Segmentação DTOS e extração de características

%WAMP (0.06)
%SSC (0.001)
%ZC (0.11)
th_ssc=0.001;
th_zc=0.11;
th_wamp = 0.06;



BL = [4.45 4.7]*1e5;

Wcrit = [0.1*fsa 1*fsa];
Wsmooth = 0.15*fsa;
k=4.0;

x = smoothdata(abs(EMGn'),'movmean',Wsmooth);
th = mean(EMGt_re(BL(1):BL(2)))+k*std(EMGt_re(BL(1):BL(2)));

[xi,xf,wT] = DTSeg(EMGn',Wcrit,x,th,T);

xi = xi(1:end-1);
xf = xf(1:end-1);
wT = wT(1:end-1);

W = [0.0 0.4]*fsa;

tic


%% 1) Preparar Matriz de EMG


%Extrair features
idx = find(ind_train==1);
for j = 1:length(idx)
 
id = idx(j);    
Sv = EMGn(xi(id):xf(id));
Sf = EMGn(xi(id)-W(1):xi(id)+W(2));

Ff(1,j)=MAV(Sf);
Ff(2,j)=WL(Sf);
Ff(3,j)=ZC(th_zc,Sf);
Ff(4,j)=SSC(th_ssc,Sf);
Ff(5,j)=WAMP(th_wamp,Sf);
Ff(6,j)=RMS(Sf);
Ff(7,j)=VAR(Sf);
Ff(8,j)=DASDV(Sf);
Ff(9,j)=IAV(Sf);
Ff(10,j)=MFL(Sf);
Ff(11,j)=MSR(Sf);
Ff(12,j)=LS(Sf,2);

Fv(1,j)=MAV(Sv);
Fv(2,j)=WL(Sv);
Fv(3,j)=ZC(th_zc,Sv);
Fv(4,j)=SSC(th_ssc,Sv);
Fv(5,j)=WAMP(th_wamp,Sv);
Fv(6,j)=RMS(Sv);
Fv(7,j)=VAR(Sv);
Fv(8,j)=DASDV(Sv);
Fv(9,j)=IAV(Sv);
Fv(10,j)=MFL(Sv);
Fv(11,j)=MSR(Sv);
Fv(12,j)=LS(Sv,2);

end


%Usar só wTv e Y_train (o test será usado a posteriori)
wTx = wT(idx);
Xtrainf = Ff;
Xtrainv = Fv;






%% 4) Classificação
%
I = eye(3);

%Avaliação do th:
%%%%%%%%%%%%%%%%
%WAMP (0.06)
%SSC (0.001)
%ZC (0.11)

%thW = (0.01:0.01:0.5);
%thW=(0:0.0005:0.005);
% thW=(0:0.01:0.2);

% for i = 1:length(thW)
%     
%     
% for j = 1:length(wT)
% 
% Sv = EMGn(xi(j):xf(j));
% Sf = EMGn(xi(j)-W(1):xi(j)+W(2));
% Fw(j)=ZC(thW(i),Sv);
% Fw(j)=ZC(thW(i),Sf);
% end
%     
% YT= wT;
% Xv = Fw';
% Xf = Fw';
% Tmat = I(YT,:);
% 
% %Não balanceado
% mdl_LDA_Tv = fitcdiscr(Xv,YT);
% evW(:,i) = kfoldLoss(crossval(mdl_LDA_Tv),'Mode','individual');
% 
% mdl_LDA_Tf = fitcdiscr(Xf,YT);
% efW(:,i) = kfoldLoss(crossval(mdl_LDA_Tf),'Mode','individual');
% 
% Xv = [];
% Xf = [];
% Fw = [];
% end
% figure
% hold on
% plot(thW,1-median(efW),' - x')
% plot(thW,1-median(evW),' - x')


%Hudgins: MAV, WL, ZC, SSC
% [1 2 3 4]
%TD4: LS, MFL, MSR, WAMP
%[12 11 10 5]
%TD9: LS, MFL, MSR, WAMP, ZC, RMS, IAV/MAV, DASDV, VAR
%[12 11 10 5 3 6 9 8 7]

YT = wTx;

for j = 1:12
    


    
Xv = Xtrainv(j,:)'; 
Xf = Xtrainf(j,:)'; 
%Tmat = I(YT,:);

%Não balanceado
mdl_LDA_Tv = fitcdiscr(Xv,YT);
cvmodelv = crossval(mdl_LDA_Tv,'CVPartition',Ck);
evi(:,j) = kfoldLoss(cvmodelv,'Mode','individual','LossFun','classiferror');

mdl_LDA_Tf = fitcdiscr(Xf,YT);
cvmodelf = crossval(mdl_LDA_Tf,'CVPartition',Ck);
efi(:,j) = kfoldLoss(cvmodelf,'Mode','individual','LossFun','classiferror');

pi(j)=ranksum(efi(:,j),evi(:,j));

Xv = [];
Xf = [];
end


% FEATURE SELECTION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fini = 1:12;
% fsel = [];
% 
% for i = 1:12
% 
% for j = 1:length(fini)
% Xv = [Fv(fsel,:);Fv(j,:)]';    
% mdl_LDA_Tv = fitcdiscr(Xv,YT);
% ev_t(:,j) = kfoldLoss(crossval(mdl_LDA_Tv),'Mode','individual');
% end
% 
% [eminv(i),inmax]=min(median(ev_t));
% fsel = [fsel fini(inmax)];
% fini(inmax) = [];
% ev_t = [];
% end
% fselv = fsel;
% 
% fini = 1:12;
% fsel = [];
% 
% for i = 1:12
% 
% for j = 1:length(fini)
% Xf = [Ff(fsel,:);Ff(j,:)]';    
% mdl_LDA_Tf = fitcdiscr(Xf,YT);
% ef_t(:,j) = kfoldLoss(crossval(mdl_LDA_Tf),'Mode','individual');
% end
% 
% [eminf(i),inmax]=min(median(ef_t));
% fsel = [fsel fini(inmax)];
% fini(inmax) = [];
% ef_t = [];
% end
% fself = fsel;
% 
% figure
% hold on
% plot(1-(eminf),'-o')
% plot(1-(eminv),'-s')
% legend('Fixed Window','Variable Length')
% 



SET{1} = [1 2 3 4];
SET{2} = [12 11 10 5];
SET{3} = [12 11 10 5 3 6 9 8 7];

for i = 1:3
%Testar com conjunto de Hudgins
%YT= wT;
Xv = Xtrainv(SET{i},:)';
Xf = Xtrainf(SET{i},:)';
%Tmat = I(YT,:);

%Não balanceado
mdl_LDA_Tv = fitcdiscr(Xv,YT);
cvmodelv = crossval(mdl_LDA_Tv,'CVPartition',Ck);
ev(:,i) = kfoldLoss(cvmodelv,'Mode','individual','LossFun','classiferror');

mdl_LDA_Tf = fitcdiscr(Xf,YT);
cvmodelf = crossval(mdl_LDA_Tf,'CVPartition',Ck);
ef(:,i) = kfoldLoss(cvmodelf,'Mode','individual','LossFun','classiferror');

Xv = [];
Xf = [];

p(i)=ranksum(ef(:,i),ev(:,i));
end








feats = {'MAV','WL','ZC','SSC','WAMP','RMS',...
     'VAR','DASDV','IAV','MFL','MSR','LS'};
 sets = {'Hudgins´','TD4','TD9'};

ei = 1-[mean(efi)' mean(evi)'];
figure
b=bar(ei*100,1);
b(1).FaceColor = [.8 .8 .8];
b(2).FaceColor = [.3 .3 1];
set(gca, 'XTickLabel',feats, 'XTick',1:numel(feats))
xtickangle(45)
legend({'Fixed Window','Variable Length'})
ylim([65 85])
ylabel('Accuracy (%)')

figure
es = 1-[mean(ef)' mean(ev)'];
b=bar(es*100,1);
b(1).FaceColor = [.8 .8 .8];
b(2).FaceColor = [.3 .3 1];
set(gca, 'XTickLabel',sets, 'XTick',1:numel(sets))
%xtickangle(45)
%legend({'Fixed Window','Variable Length'})
ylim([65 90])
ylabel('Accuracy (%)')

% %% Plotar resultados
% %
% a = 0.25; %distância entre box
% bW = 0.15; %largura do box
% position = [];
% tickspos=[];
% color = ['y', 'm', 'c']; %cor do box
% C=[];
% for i = 1:J
%    position = cat(2,position,[i-a i i+a]); 
%    tickspos = cat(2,tickspos,mean(position(1+(i-1)*K:3+(i-1)*K))); 
%    C = cat(2,C,color);
% end
% 
% figure
% boxplot([(1-e(:,:,1)) (1-e(:,:,2)) (1-e(:,:,3))],...
%     [ K*(1:J)-2 K*(1:J)-1 K*(1:J)],'widths',bW,'positions',position,'Colors','kkk')
% 
% set(gca,'xtick',tickspos)
% set(gca,'xticklabel',{'600 ms','800 ms','1 s','2 s','3 s','5s','10s'})
% 
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%    patch(get(h(j),'XData'),get(h(j),'YData'),C(j),'FaceAlpha',.5);
% end
% 
% c = get(gca, 'Children');
% 
% hleg1 = legend(c(1:3), '5 \sigma','10 \sigma','15 \sigma');
% 
% set(gca,'xtick',tickspos)
% set(gca,'xticklabel',{'600 ms','800 ms','1 s','2 s','3 s','5s','10s'})
% 


%% Funções
% 




function [xi,xf,wT] = DTSeg(data,Wcrit,x,th,T)
%Double Threshold segmentation (variable length)
% A.K.A Onset Segmentation
%
% -> Find segments which are:          
%    *Longer than Wcrit               
%    *Greater than th                  
%
% Inputs:
% - data: data vector
% - win: Signal smooth window size
% - Wcrit: Time threshold (min segmentation length)
% - th: Amplitude threshold
% - method: smoothing method
% Outputs:
% - segM: segment cell
% - where: segment extension

%data length
N= max(size(data));
%smooth absolute signal

%PB = fH*2/fs;
%[Bb,Ab] = butter(6,PB,'low');
%x=filter(Bb,Ab,abs(ex));


%th = mean(x(BL(1):BL(2)))+k*std(x(BL(1):BL(2)));

%th = k*std(x(BL(1):BL(2)));

%segM=[];
%where=zeros(1,N);
wT = [];

i=1;
j=1;
while i~=N

   %Signal is under th... 
   while abs(x(i)<th) && i~=N 
       i=i+1;        
   end
   %...signal is over th!
   
   %Since Ci, there is a new segment...
   Ci=i;
   
   %...still over th....
   while abs(x(i)>th) && i~=N  
       i=i+1;       
   end
   %...and its gone. Now it is under th.
    
   %The segment ends at Cf. 
   Cf=i;
   
   %The segment existance is limited by Ci and Cf.
   %The difference is the segment length 
   C=Cf-Ci;
   
   %Is this segment longer than Wcrit?
   if (C>Wcrit(1)) && (C<Wcrit(2))
       
       xi(j) = Ci;
       xf(j) = Cf;
       %We got a new segment!
       w=Ci:Cf;             %segment boundary
       %w_fix=Ci:(Ci+W-1)
       
       Tw = T(w,:);
       [~,wi]=max(sum(Tw));
       wT=cat(2,wT,wi);
       
       %where(w) = 1;        %If there is a segment, there is a 1.
       %segM(j,:) = data(w);   %Add the segment to the cell
       %segM = data(w);   %Add the segment to the cell
       
       
       
       j=j+1;
       
       
    
   end
   
end

end

