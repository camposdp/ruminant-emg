% Daniel Prado de Campos - UTFPR - 22/04/2019
% Avaliação estatística - DTOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 0) Carregamento

clear
clc
close all


load EMG_filt_e_norm
load EMG_visual
%load EMG_segment_labels
load cont
load Error_seg_eval

fsa = 2000;

Xi = xi;
Xf = xf;
X = (Xf - Xi)/fsa*1000;
EMGt_re = EMGt(1:xf(end));





Xi = xi;
Xf = xf;
X = (Xf - Xi)/fsa*1000;

%% 1) Plotar dados da avaliação do método (confusao_segment_visual_DTOS)

figure
plot(Kvar,100*c,'x-','LineWidth',0.8);
xlabel('k (threshold = \mu + k*\sigma)')
ylabel('Error (%)')
legendCell = cellstr(num2str(1000*Wvar', 'W=%d ms'));
legend(legendCell,'Location','SouthEast')
legend('boxoff')
%xticks(Kvar)

[cmin,argmin] = min(c);
axes('Position',[.58 .68 .3 .2])
box on
plot(1000*Wvar,cmin*100,'-or','LineWidth',1)
title('Minimum error')
xlabel('Smooth Window Size (ms)')
ylabel('Error (%)')
xlim([0 0.17]*1000)
ylim([7 9])
xticks(1000*Wvar)
xtickangle(45)
tickCell = cellstr(num2str(Kvar(argmin)', 'k=%1.1f'));
xt = Wvar*1000;
xt(3:end) = xt(3:end)+3;
yt = cmin*100;
yt(3:end) = yt(3:end)-0.2;
yt(1:2) = yt(1:2)+0.2;
text(xt',yt',tickCell,'FontSize',8);


%% 2) Análise visual

% Xi = [];
% Xf = [];
% for i = 1:13
%     s = ['cont' num2str(i)];
%     load(s)
%     Xi = [Xi xi];
%     Xf = [Xf xf];
%     
%     clear xf
%     clear xi
% end




cov(1,1,1)=5;
mu = 300;
St = struct('mu',mu,'Sigma',cov);

cov(1,1,1:2)=[5 5];
mu = [100 300]';
St2 = struct('mu',mu,'Sigma',cov);

figure
subplot(1,2,1)
[Px]=getGausMix1(X,St,true);
xlim([0 500])
ylim([0 6e-3])
title('Visual Segmentation')
%% 4) Segmentação automática

BL = [4.45 4.7]*1e5;
Wcrit = [0.1*fsa 1*fsa]; 
%Wcrit = [0.015*fsa 0.6*fsa];
W = 0.03*fsa;
k=2.5;
%W=0.05*fsa;

%[~,ex]=energyop(EMGt_re,false);
%EMGs = smooth(abs(ex),W,'lowess');
x = smoothdata(abs(EMGt_re),'movmean',W);

th = mean(x(BL(1):BL(2)))+k*std(x(BL(1):BL(2)));

% figure
% hold all
% plot(EMGt_re)
% plot(x)
% plot(th*ones(length(x),1))

[Yi,Yf] = DTSeg(EMGt_re, Wcrit,x,th);
Y = (Yf - Yi)/fsa*1000;

subplot(1,2,2)
[Py]=getGausMix(Y,St2,true);
xlim([0 500])
ylim([0 6e-3])  
title('Automatic Segmentation (DTOS)')

Xbin = zeros(length(EMGt_re),1);
Ybin = zeros(length(EMGt_re),1);

for i = 1:length(X)
    Xbin(round(Xi(i)):round(Xf(i)))=1;
end

for i = 1:length(Y)
    Ybin(round(Yi(i)):round(Yf(i)))=1;
end
%Xbin = zeros(length(EMGt_re),1);
[c,~,~,per] = confusion(Xbin(1:end-1)',Ybin');


figure
subplot(1,2,1)
hold all
xi=(0:1e-4:600);   

%nx2 = normpdf(xi,Px(2),(Px(4))^0.5)*Px(6);
nx1 = normpdf(xi,Px(1),(Px(2))^0.5)*Px(3);

%Q = nx1 + nx2;
%Q2 = nx2;   
Q = nx1;

%plot(xi,nx1,'r-','LineWidth',1)
plot(xi,Q,'k--','LineWidth',1)


Wv = Wvar;
kv =Kvar(argmin);
Nw = length(Wv);
for i = 1:Nw
W = Wv(i)*fsa;
k = kv(i);
x = smoothdata(abs(EMGt_re),'movmean',W);
th = mean(x(BL(1):BL(2)))+k*std(x(BL(1):BL(2)));    
[Yi,Yf] = DTSeg(EMGt_re, Wcrit,x,th);
 Y = (Yf - Yi)/fsa*1000;    
[Py]=getGausMix(Y,St2,false);

Ybin = zeros(length(EMGt_re),1);
for q = 1:length(Y)
    Ybin(round(Yi(q)):round(Yf(q)))=1;
end

[c2(i),~,~,per2] = confusion(Xbin(1:end-1)',Ybin');
perV(i,:)=per2(1,:);


ny1 = normpdf(xi,Py(1),(Py(3))^0.5)*Py(5);
ny2 = normpdf(xi,Py(2),(Py(4))^0.5)*Py(6);



P(:,i) = ny2+ny1;
P2(:,i) = ny2;

u(i) = Py(1);

%plot(xi,ny1,'-.')
h=plot(xi,ny2+ny1,'-','Linewidth',0.8);
hc(i,:)=h.Color;
end

du = abs(u-Px(1))';

legendCell = cellstr([num2str(1000*Wv', 'W=%d ms') num2str(kv','\t,k=%1.1f')]);
legend(['Visual Reference';legendCell])
ylim([0 8e-3]);
xlabel('sEMG Period (ms)')
ylabel('Probability Density')
lgd=legend('boxoff')
lgd.FontSize = 8;

ms = 100;
%axes('Position',[.20 .68 .2 .22])
%box on
% subplot(2,2,3)
% hold on
% title('Main Component PDF Mean')
% scatter(Wv*1000,u,ms,hc,'filled')
% plot(Wv*1000,Px(1)*ones(Nw,1),'k--')
% %ylim([150 350])
% %xlim([0 160])
% xticks(Wv*1000)
% xtickangle(45)
% xlabel('Smoothing Window (ms)')
% ylabel('Mean (ms)')
% %legend({'Automatic Segmentation','Reference'},...
% %    'Location','SouthEast')
% %legend('boxoff')
% grid on

dist=KLDiv(P',Q);
%dist2=KLDiv(P2',Q2);

subplot(1,2,2)
hold all
scatter(Wv*1000,dist,ms,hc,'fill','^')
%scatter(Wv*1000,dist2,ms,hc,'fill','s')
xticks(Wv*1000)
xlim([0 160])
ylim([0 0.25])
xticks(Wv*1000)
xtickangle(45)
xlabel('Smoothing Window (ms)')
ylabel('Kullback-Leibler divergence')
%legend({'Mixture','Main distribution'})
title('GMM KLD')
    grid on
    
% subplot(2,3,6)
% hold all
% %scatter(Wv*1000,dist,ms,hc,'fill','^')
% scatter(Wv*1000,dist2,ms,hc,'fill','v')
% xticks(Wv*1000)
% xlim([0 160])
% xticks(Wv*1000)
% xtickangle(45)
% xlabel('Smoothing Window (ms)')
% ylabel('Kullback-Leibler divergence')
% %legend({'Mixture','Main distribution'})
% title('Main Component PDF KLD')
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions
%

function [xi,xf] = DTSeg(data,Wcrit,x,th)
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
%wT = [];

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
       %w=Ci:(Ci+W-1);             %segment boundary
       
       %Tw = T(w,:);
       %[~,wi]=max(sum(Tw));
       %wT=cat(2,wT,wi);
       
       %where(w) = 1;        %If there is a segment, there is a 1.
       %segM(j,:) = data(w);   %Add the segment to the cell
       %segM = data(w);   %Add the segment to the cell
       
       
       
       j=j+1;
       
       
    
   end
   
end

end

function [P]=getGausMix1(X,St,flag)





opt = struct('MaxIter',3000);
GMModel = fitgmdist(X',1,'Start',St,'Options',opt);


 
u1 = GMModel.mu(1);
%u2 = GMModel.mu(2);
sigma1 = GMModel.Sigma(1);
%sigma2 = GMModel.Sigma(2);
p1 = GMModel.ComponentProportion(1);
%p2 = GMModel.ComponentProportion(2);

%P = [u1 u2 sigma1 sigma2 p1 p2];
P = [u1 sigma1 p1];

if flag == logical(1)

%figure
hold all


hGM = histogram(X,'Normalization','pdf','BinWidth',25)
hGM.FaceColor = [0.8 0.8 0.8]
hGM.EdgeColor = [1 1 1]

eGM = ezplot(@(x)pdf(GMModel,[x]),[0 max(X)*0.9]);
eGM.LineWidth=1.0;
eGM.Color=[0 0 0];

x=(0:1e-4:max(X));   
y1 = normpdf(x,u1,(sigma1)^0.5)*p1;
%y2 = normpdf(x,u2,(sigma2)^0.5)*p2;
%plot(x,y1,'r--','LineWidth',1)
%plot(x,y2,'b--','LineWidth',1)
ylabel('Probability Density')
xlabel('Period (ms)')
s1 = ['\mu = ' num2str(u1,'%10.1f')...
      ', \sigma = ' num2str(sigma1^0.5,'%10.1f')...
      ];
%s2 = ['\mu_{2} = ' num2str(u2,'%10.1f')...
%      ', \sigma_{2} = ' num2str(sigma2^0.5,'%10.1f')...
%      ', p_{2} = ',num2str(p2)];  
%legend({'Histogram','Gaussian Mixture Fit',s1,s2},'FontSize',8)
%title('sEMG Period Probability Density Function Fit')
title('')
end

end


function [P]=getGausMix(X,St,flag)





opt = struct('MaxIter',3000);
GMModel = fitgmdist(X',2,'Start',St,'Options',opt);



u1 = GMModel.mu(1);
u2 = GMModel.mu(2);
sigma1 = GMModel.Sigma(1);
sigma2 = GMModel.Sigma(2);
p1 = GMModel.ComponentProportion(1);
p2 = GMModel.ComponentProportion(2);

P = [u1 u2 sigma1 sigma2 p1 p2];


if flag == logical(1)

%figure
hold all


hGM = histogram(X,'Normalization','pdf','BinWidth',25)
hGM.FaceColor = [0.8 0.8 0.8]
hGM.EdgeColor = [1 1 1]

x=(0:1e-4:max(X));   
y1 = normpdf(x,u1,(sigma1)^0.5)*p1;
y2 = normpdf(x,u2,(sigma2)^0.5)*p2;

eGM = ezplot(@(x)pdf(GMModel,[x]),[0 max(X)*0.9]);
eGM.LineWidth=1.0;
eGM.Color=[0 0 0];


plot(x,y1,'r--','LineWidth',1)
plot(x,y2,'b--','LineWidth',1)
ylabel('Probability Density')
xlabel('Period (ms)')
s1 = ['\mu_{1} = ' num2str(u1,'%10.1f')...
      ', \sigma_{1} = ' num2str(sigma1^0.5,'%10.1f')...
      ', p_{1} = ',num2str(p1)];
s2 = ['\mu_{2} = ' num2str(u2,'%10.1f')...
      ', \sigma_{2} = ' num2str(sigma2^0.5,'%10.1f')...
      ', p_{2} = ',num2str(p2)];  
%legend({'Histogram','Gaussian Mixture Fit',s1,s2},'FontSize',8)
%title('sEMG Period Probability Density Function Fit')
title('')
end

end



function dist=KLDiv(P,Q)
%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1
if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end
if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!') 
end
% normalizing the P and Q
if size(Q,1)==1
    Q = Q ./sum(Q);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./repmat(Q,[size(P,1) 1]));
    temp(isnan(temp))=0;% resolving the case when P(i)==0
    dist = sum(temp,2);
    
    
elseif size(Q,1)==size(P,1)
    
    Q = Q ./repmat(sum(Q,2),[1 size(Q,2)]);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./Q);
    temp(isnan(temp))=0; % resolving the case when P(i)==0
    dist = sum(temp,2);
end

end
