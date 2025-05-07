%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daniel Prado de Campos 24/05/2019
% FDDL Wrapper
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
% Ck = cvpartition(Ninst,'KFold',10);
rng('default');

[~,Ns] = size(SS);
[~,Nr] = size(SR);


% randIdxSS = randi(Ns,Ninst,1);
% randIdxSR = randi(Nr,Ninst,1);
SRrand = SS(:,randIdxSS);
SSrand = SR(:,randIdxSR);

label_train = [ones(900,1);2*ones(900,1)];
label_test = [ones(100,1);2*ones(100,1)];

%% FDDL
fprintf('Step 1.0. Iniciando FDDL...\n');

lambda1 = 0.01;
lambda2 = 0.01;
Kdic = 400;


SNR = (0:20);

for i = 1:length(SNR)


SRnoisy = awgn(SRrand,SNR(i),'measured');
SSnoisy = awgn(SSrand,SNR(i),'measured');

for j = 1:10
  

trIdx = Ck.training(j);
teIdx = Ck.test(j);

%Treinar com ruído
%Train1 = SRnoisy(:,trIdx );
%Train2 = SSnoisy(:,trIdx );

%Treinar sem ruído
Train1 = SRrand(:,trIdx );
Train2 = SSrand(:,trIdx );

Test1 = SRnoisy(:,teIdx);
Test2 = SSnoisy(:,teIdx);   

Y_train_k = [Train1 Train2];
Y_test_k = [Test1 Test2];

[AC, rt, Di, CoefMi, optsi] = myFDDL(Y_train_k, label_train', Y_test_k, label_test', ...
                            Kdic, lambda1, lambda2);

acc(i,j) = max(AC)
D{i,j} = Di;
CoefM{i,j} = CoefMi;
opts{i,j}=optsi;

Y_train_k = [];  
Y_test_k = [];    

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



function [AC, rt, D, CoefM, opts] = myFDDL(Y_train, label_train, Y_test , label_test, ...
                            k, lambda1, lambda2)
% function [acc, rt] = FDDL_wrapper(Y_train, label_train, Y_test , ...
%       label_test, k, lambda1, lambda2)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    if nargin == 0 % test mode
        dataset = 'myYaleB';
        N_train = 10;        
        [~, Y_train, Y_test, label_train, label_test] = ...
            train_test_split(dataset, N_train);        
        k = 8;
        lambda = 0.001;
        eta = 0.01;
    end 
    C                = max(label_train);
    k0               = 0;    
    opts.k           = k;
    opts.k0          = 0;
    opts.show_cost   = 0;
    opts.lambda1     = lambda1;
    opts.lambda2     = lambda2;
    opts.lambda3     = 0;
    opts.D_range     = k*(0:C);
    opts.D_range_ext = [opts.D_range k*C+k0];
    opts.initmode    = 'normal';   
    opts.max_iter    = 100;
    opts             = initOpts(opts);
    opts.verbose      = true;
    opts.tol         = 1e-8;
    %% Train 
    fprintf('Step 1.1. Treinando dicionário...\n');
    [D, ~, ~, ~, CoefM, ~, opts, rt] = ...
                    LRSDL(Y_train, label_train, opts);
    Y_range = label_to_range(label_train);
    C = max(label_train);
    opts.verbose = 0;
    opts.weight = 0.1;
    acc = [];
   % for vgamma = [0.0001, 0.001, 0.01, 0.1]
    fprintf('Step 1.2. Encontrando gama...\n');
    vg = [(0.00005:0.00005:0.00045) (0.0005:0.0005:0.01) (0.02:0.01:1)];
    for vgamma = vg
        opts.gamma = vgamma;
        pred = FDDL_pred(Y_test, D, CoefM, opts);
        acc1 = double(numel(find(pred == label_test)))/...
            numel(label_test);
        fprintf('gamma = %.4f, acc = %.4f\n', vgamma, acc1);
        acc = [acc acc1];
    end 
    [AC,qual] = max(acc);
    
    
    
    best_gamma = vg(qual);

    opts.gamma = best_gamma;

 
    
    
end 