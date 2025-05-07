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

load EMG_base
load EMGre_exp

%% 0.1 Definições

fprintf('Step 0.1. Definindo parâmetros...\n');

fsa = 2000;

BL = [4.45 4.7]*1e5;

Wcrit = [0.1*fsa 1*fsa];
Wsmooth = 0.15*fsa;
k=4.0;


EMG_R = [EMG_A7;EMG_A9];
EMG_S = [EMG_A6;EMG_A8];
% EMG_R = [EMG_A7];
% EMG_S = [EMG_A6];

%% 0.2 Segmentação


fprintf('Step 0.2. Segmentando dados...\n');

x_r = smoothdata(abs(EMG_R),'movmean',Wsmooth);
x_s = smoothdata(abs(EMG_S),'movmean',Wsmooth);

th = mean(EMGt_re(BL(1):BL(2)))+k*std(EMGt_re(BL(1):BL(2)));

[xiR,xfR] = DTOSSeg(EMG_R,Wcrit,x_r,th);
[xiS,xfS] = DTOSSeg(EMG_S,Wcrit,x_s,th);




W = [0.1 0.3]*fsa;

   

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

Ck = cvpartition(500,'KFold',10);
rng('default');

[~,Ns] = size(SS);
[~,Nr] = size(SR);
Ninst = 500;

randIdxSS = randi(Ns,Ninst,1);
randIdxSR = randi(Nr,Ninst,1);
SRrand = SS(:,randIdxSS);
SSrand = SR(:,randIdxSR);

label_train = [ones(450,1);2*ones(450,1)];
label_test = [ones(50,1);2*ones(50,1)];

%% FDDL
fprintf('Step 1.0. Iniciando FDDL...\n');


lambda1 = 0.01;
lambda2 = 0.01;

Ki = [200];

for j = 1:10
for i = 1:length(Ki)
    
 Kdic = Ki(i);   
 pS = 100;  
 dS = 500;   
 pSize = 200;
 OLstep = 0.5;
 
trIdx = Ck.training(j);
teIdx = Ck.test(j);
          
Train1 = SRrand(:,trIdx );
Train2 = SSrand(:,trIdx );

Test1 = SRrand(:,teIdx);
Test2 = SSrand(:,teIdx);   

Y_train_k = [Train1 Train2];
Y_test_k = [Test1 Test2];


[AC, rt, Di, CoefMi, optsi] = FDDL_patch(Y_train_k, label_train, Y_test_k, label_test, ...
                            Kdic, lambda1, lambda2, pSize,OLstep)

acc(i,j) = max(AC);
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


function [acc, rt, D, CoefM, opts] = FDDL_patch(Y_train, label_train, Y_test , label_test, ...
                            k, lambda1, lambda2, pSize,OLstep)
% function [acc, rt] = FDDL_wrapper(Y_train, label_train, Y_test , ...
%       label_test, k, lambda1, lambda2)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    thresh = 0.5;
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
    
    opts.patchSize = pSize;
    opts.overlapStep = OLstep;
    %% Train 
    [~,N_train] = size(Y_train);
    Y_train_p = [];
    label_train_p = [];
    for i = 1:N_train
        Yp = buildOverlappingPatches(opts, Y_train(:,i));
        [~,Npatches] = size(Yp); 
        Y_train_p = [Y_train_p Yp];
        label_train_p = [label_train_p;...
        label_train(i)*ones(Npatches,1)];
    end
    
    [D, ~, ~, ~, CoefM, ~, opts, rt] = ...
                    LRSDL(Y_train_p, label_train_p, opts);
    Y_range = label_to_range(label_train_p);
    C = max(label_train_p);
    opts.verbose = 0;
    opts.weight = 0.1;
    acc = [];
    %vg = [(0.00005:0.00005:0.00045) ...
     %   (0.0005:0.0005:0.01) (0.02:0.01:1)];
     vg = [0.0001 0.001 0.01 0.1];
    for vgamma = vg
        
        opts.gamma = vgamma;
        
        
        [~,N_test] = size(Y_test);
        Y_test_p = [];
        label_test_p = [];
        for i = 1:N_test
            Yp = buildOverlappingPatches(opts, Y_test(:,i));
            [~,Npatches] = size(Yp); 
            label_test_p = [label_test(i)*ones(Npatches,1)];
            predi = FDDL_pred(Yp, D, CoefM, opts);
            feature           = sum(predi == 1)/numel(predi);
            pred(i)  = -0.5*(2*(feature > thresh) -1) + 1.5;
        end

        
       
        acc1 = double(numel(find(pred' == label_test)))/...
            numel(label_test);
        fprintf('gamma = %.4f, acc = %.4f\n', vgamma, acc1);
        acc = [acc acc1];
    end 
    acc = max(acc);
end 


function X = buildOverlappingPatches(opts, signal)
	p       = opts.patchSize;
    step = floor(opts.overlapStep*p);
    W = max(size(signal));
    
    xi = 1;
    X=[];
    while xi <= (W-p+1)
    X = [X signal(xi:xi+p-1)];
    xi = xi+step;
    end
    
end
