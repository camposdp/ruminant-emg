%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daniel P. Campos - UTFPR/CPGEI - 24/05/2019
% Segmentar, extrair features e classificar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Abrir Dados
clear
close all
clc


load EMG_base_exp
load EMGre_exp
load KFOLD_1000 

%load KFOLD_500

%% Definições

fsa = 2000;

BL = [4.45 4.7]*1e5;

Wcrit = [0.1*fsa 1*fsa];
Wsmooth = 0.03*fsa;
k=2.5;


EMG_R = [EMG_A7;EMG_A9];
EMG_S = [EMG_A6;EMG_A8;EMG_A10;EMG_A11;EMG_A12;EMG_A13;EMG_A14];

%% Segmentação

x_r = smoothdata(abs(EMG_R),'movmean',Wsmooth);
x_s = smoothdata(abs(EMG_S),'movmean',Wsmooth);

th = mean(EMGt_re(BL(1):BL(2)))+k*std(EMGt_re(BL(1):BL(2)));

[xiR,xfR] = DTOSSeg(EMG_R,Wcrit,x_r,th);
[xiS,xfS] = DTOSSeg(EMG_S,Wcrit,x_s,th);


%% Extração de features

W = [0.1 0.4]*fsa;
%SR = [];

SNR = (0:20);

for isnr = 1:length(SNR)


Nr = length(xiR);
for id = 1:Nr
 
SvR = EMG_R(xiR(id):xfR(id));
SfR = EMG_R(xiR(id)-W(1):xiR(id)+W(2)-1);

%SR = [SR SfR];

[FvR(:,id)] = getFeaturesNoisy(SvR,SNR(isnr));
[FfR(:,id)] = getFeaturesNoisy(SfR,SNR(isnr));

[FvR0(:,id)] = getFeaturesClean(SvR);
[FfR0(:,id)] = getFeaturesClean(SfR);



end



Ffr1(isnr,:) = FfR(1,:);
Ffr2(isnr,:) = FfR(2,:);
Ffr3(isnr,:) = FfR(3,:);
Ffr4(isnr,:) = FfR(4,:);


%SS = [];
Ns = length(xiS);
for id = 1:Ns

SvS = EMG_S(xiS(id):xfS(id));
SfS = EMG_S(xiS(id)-W(1):xiS(id)+W(2)-1);

[FvS(:,id)] = getFeaturesNoisy(SvS,SNR(isnr));
[FfS(:,id)] = getFeaturesNoisy(SfS,SNR(isnr));

[FvS0(:,id)] = getFeaturesClean(SvS);
[FfS0(:,id)] = getFeaturesClean(SfS);

%SS = [SS SfS];

end

Ffs1(isnr,:) = FfS(1,:);
Ffs2(isnr,:) = FfS(2,:);
Ffs3(isnr,:) = FfS(3,:);
Ffs4(isnr,:) = FfS(4,:);



%% Balanceamento 
rng('default');
Ninst = 1000;


%nrb = randi(Nr,Ninst,1);
%nsb = randi(Ns,Ninst,1);
YT = [ones(Ninst,1);2*ones(Ninst,1)];
Fv = [FvR(:,randIdxSR) FvS(:,randIdxSS)];
Ff = [FfR(:,randIdxSR) FfS(:,randIdxSS)];

Fv0 = [FvR0(:,randIdxSR) FvS0(:,randIdxSS)];
Ff0 = [FfR0(:,randIdxSR) FfS0(:,randIdxSS)];

%YT = [ones(Nr,1);2*ones(Ns,1)];
% Fv = [FvR FvS];
% Ff = [FfR FfS];


for i = 1:10
Ck2test(:,i) = [Ck.test(i);Ck.test(i)];
Ck2training(:,i) = [Ck.training(i);Ck.training(i)];
end


%% Definição de Sets



SET{1} = [1 2 3 4];
SET{2} = [1 2 3 4 6 13 14 15 16 17 18];
SET{3} = [12 11 10 5];
SET{4} = [12 11 10 5 3 6 9 8 7];


%% Classificação


for set = 1:length(SET)
%Testar com conjunto de Hudgins
%YT= wT;
Xv = Fv(SET{set},:)';
Xf = Ff(SET{set},:)';



Xv0 = Fv0(SET{set},:)';
Xf0 = Ff0(SET{set},:)';

% Validação cruzada
% 
% mdl_LDA_Tv = fitcdiscr(Xv,YT);
% cvmodelv = crossval(mdl_LDA_Tv,'partition',Ck);
% ev(:,i) = kfoldLoss(cvmodelv,'Mode','individual','LossFun','classiferror');
% 
% mdl_LDA_Tf = fitcdiscr(Xf,YT);
% cvmodelf = crossval(mdl_LDA_Tf,'partition',Ck);
% ef(:,i) = kfoldLoss(cvmodelf,'Mode','individual','LossFun','classiferror');

    for j = 1:10

    idxTr = Ck2training(:,j);
    idxTe = Ck2test(:,j);

    mdl_LDA_Tv = fitcdiscr(Xv(idxTr,:),YT(idxTr));
    ev = loss(mdl_LDA_Tv,Xv0(idxTe,:),YT(idxTe),...
        'LossFun','classiferror');

    mdl_LDA_Tf = fitcdiscr(Xf(idxTr,:),YT(idxTr));
    ef = loss(mdl_LDA_Tf,Xf0(idxTe,:),YT(idxTe),...
        'LossFun','classiferror');

    if set == 1
        ev1(isnr,j) = ev; 
        ef1(isnr,j) = ef; 
    elseif set ==2
        ev2(isnr,j) = ev; 
        ef2(isnr,j) = ef; 
    elseif set ==3
        ev3(isnr,j) = ev; 
        ef3(isnr,j) = ef; 
    elseif set ==4
        ev4(isnr,j) = ev; 
        ef4(isnr,j) = ef; 
    end



    end

Xv = [];
Xf = [];

%p(i)=ranksum(ef(:,i),ev(:,i));
end
end

plot(SNR,[mean(1-ef1');mean(1-ef2');mean(1-ef3');mean(1-ef4')]','-x')
xticks(SNR)
legend('MS1','MS2','TD4','TD9')

%boxplot([1-ev 1-ef])


%% Plotar 


% feats = {'MAV','WL','ZC','SSC','WAMP','RMS',...
%      'VAR','DASDV','IAV','MFL','MSR','LS'};
%  sets = {'Hudgins´','TD4','TD9'};
% 
%  pcrit = 0.05;
%  
% ei = 1-[mean(efi)' mean(evi)'];
% figure
% hold all
% b1=bar(ei*100,1);
% 
% b1(1).FaceColor = [.8 .8 .8];
% b1(2).FaceColor = [.3 .3 1];
% set(gca, 'XTickLabel',feats, 'XTick',1:numel(feats))
% 
% xtickangle(45)
% 
% 
% ylim([50 85])
% ylabel('Accuracy (%)')
% getSigstar(b1,ei,pi,pcrit)
% %legend({'Fixed Window','Variable Length'},'Location','NorthWest')
% 
% figure
% hold all
% es = 1-[mean(ef)' mean(ev)'];
% b2=bar(es*100,1);
% 
% b2(1).FaceColor = [.8 .8 .8];
% b2(2).FaceColor = [.3 .3 1];
% set(gca, 'XTickLabel',sets, 'XTick',1:numel(sets))
% %xtickangle(45)
% 
% 
% ylim([65 85])
% ylabel('Accuracy (%)')
% getSigstar(b2,es,p,pcrit)
% legend({'Fixed Window','Variable Length'},'Location','NorthWest')
%% Funnções



function [F] = getFeaturesClean(S)



th_ssc=0.001;
th_zc=0.11;
th_wamp = 0.06;

F(1)=MAV(S);
F(2)=WL(S);
F(3)=ZC(th_zc,S);
F(4)=SSC(th_ssc,S);
F(5)=WAMP(th_wamp,S);
F(6)=RMS(S);
F(7)=VAR(S);
F(8)=DASDV(S);
F(9)=IAV(S);
F(10)=MFL(S);
F(11)=MSR(S);
F(12)=LS(S,2);
ar = AR6(S);
F(13)=ar(1);
F(14)=ar(2);
F(15)=ar(3);
F(16)=ar(4);
F(17)=ar(5);
F(18)=ar(6);


end


function [F] = getFeaturesNoisy(Sclean,SNR)

S = awgn(Sclean,SNR,'measured');


th_ssc=0.001;
th_zc=0.11;
th_wamp = 0.06;

F(1)=MAV(S);
F(2)=WL(S);
F(3)=ZC(th_zc,S);
F(4)=SSC(th_ssc,S);
F(5)=WAMP(th_wamp,S);
F(6)=RMS(S);
F(7)=VAR(S);
F(8)=DASDV(S);
F(9)=IAV(S);
F(10)=MFL(S);
F(11)=MSR(S);
F(12)=LS(S,2);
ar = AR6(S);
F(13)=ar(1);
F(14)=ar(2);
F(15)=ar(3);
F(16)=ar(4);
F(17)=ar(5);
F(18)=ar(6);


end


function [xi,xf] = DTOSSeg(data,Wcrit,x,th)
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
%       w=Ci:Cf;             %segment boundary
       %w_fix=Ci:(Ci+W-1)
       
%        Tw = T(w,:);
%        [~,wi]=max(sum(Tw));
%        wT=cat(2,wT,wi);
       
       %where(w) = 1;        %If there is a segment, there is a 1.
       %segM(j,:) = data(w);   %Add the segment to the cell
       %segM = data(w);   %Add the segment to the cell
       
       
       
       j=j+1;
       
       
    
   end
   
end

end