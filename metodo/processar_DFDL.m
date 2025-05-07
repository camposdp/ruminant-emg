%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daniel Prado de Campos 24/05/2019
% DFDL Wrapper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0.0 Abrir Dados
clear
close all
clc
tic
fprintf('Step 0.0. Abrindo arquivos...\n');

load EMG_base_exp
load EMGre_exp
load KFOLD_1000

%% 0.1 Definições

fprintf('Step 0.1. Definindo parâmetros...\n');

fsa = 2000;

BL = [4.45 4.7]*1e5;

Wcrit = [0.1*fsa 1*fsa];
Wsmooth = 0.03*fsa;
k=2.5;


EMG_R = [EMG_A7;EMG_A9];
EMG_S = [EMG_A6;EMG_A8;EMG_A10;EMG_A11;EMG_A12;EMG_A13;EMG_A14];
% EMG_R = [EMG_A7];
% EMG_S = [EMG_A6];

%% 0.2 Segmentação


fprintf('Step 0.2. Segmentando dados...\n');

x_r = smoothdata(abs(EMG_R),'movmean',Wsmooth);
x_s = smoothdata(abs(EMG_S),'movmean',Wsmooth);

th = mean(EMGt_re(BL(1):BL(2)))+k*std(EMGt_re(BL(1):BL(2)));

[xiR,xfR] = DTOSSeg(EMG_R,Wcrit,x_r,th);
[xiS,xfS] = DTOSSeg(EMG_S,Wcrit,x_s,th);




W = [0.1 0.4]*fsa;

   

SR = [];
Nr = length(xiR);
for id = 1:Nr-1
 
SvR = EMG_R(xiR(id):xfR(id));
SfR = EMG_R(xiR(id)-W(1):xiR(id)+W(2)-1);

SR = [SR SfR];

end

SS = [];
Ns = length(xiS);
for id = 1:Ns-1

SvS = EMG_S(xiS(id):xfS(id));
SfS = EMG_S(xiS(id)-W(1):xiS(id)+W(2)-1);


SS = [SS SfS];

end


%% 0.3 Balanceamento 


fprintf('Step 0.3. Balanceando instâncias...\n');
Ninst = 1000;
%Ck = cvpartition(Ninst,'KFold',10);
rng('default');

[~,Ns] = size(SS);
[~,Nr] = size(SR);


%randIdxSS = randi(Ns,Ninst,1);
%randIdxSR = randi(Nr,Ninst,1);
SRrand = SS(:,randIdxSS);
SSrand = SR(:,randIdxSR);

label_train = [ones(900,1);2*ones(900,1)];
label_test = [ones(100,1);2*ones(100,1)];

%% DFDL

 
 
 Ki = [50 100 200 300 400 500];

for i = 1:length(Ki)
for j = 1:10
    
 pS = 100;  
 dS = Ki(i); 
 oS = 0.5;
 
trIdx = Ck.training(j);
teIdx = Ck.test(j);
          
Train1 = SRrand(:,trIdx );
Train2 = SSrand(:,trIdx );

Test1 = SRrand(:,teIdx);
Test2 = SSrand(:,teIdx);   
    

% Parameters - Patches
pars.overlapStep = oS; 
pars.numPatches = 15000; 
pars.patchSize = pS;

% Parameters DFDL
pars.dictsize = dS; %Dicionary size
pars.K        = pars.dictsize*ones(1,2); 
                      % pars.K: #bases/class
pars.rho = 0.1;                      
pars.max_iter = 50;	  % 
pars.lambda   = 0.1;  % , see TMI paper
pars.gamma    = 0.1;  % , see TMI paper

% Parameters - OMP
paramOMP.eps        = 1e-5;
paramOMP.numThreads = -1;  
pars.paramOMP       = paramOMP;

% Não plotar ROC
roc = false;

acc(i,j) = DFDLwrapper(Train1,Train2,Test1,Test2,pars,roc)

Train1 = [];
Train2 = [];
Test1 = [];
Test2 = [];


t(i,j)=toc
end
end


%% Funções

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