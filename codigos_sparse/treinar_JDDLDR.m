clear all;clc;
addpath('utilities');



clear
clc
close all

%load Dic_kSVD
load EMG_labels_and_windows
load labels_val_e_test


fs = 2000;

Nb = 0.2*fs;
Nf = 0.6*fs;
w0 = Nb+Nf;

ttls = wTt;
trls = wTv;


N_t = length(ttls);
N_v = length(trls);


%% 1) Preparar Matriz de EMG

d = 1;
w = w0/d;
Train_DAT = zeros(w,N_v);
Test_DAT = zeros(w,N_t);

for i = 1:N_v
Train_DAT(:,i) = decimate(EMGn(xiv(i)-Nb:xiv(i)+Nf-1),d); 
end

for i = 1:N_t
Test_DAT(:,i) = decimate(EMGn(xit(i)-Nb:xit(i)+Nf-1),d); 
end


numcomps = 500;      % tr_num = 3;
nDim     = 4;    % feature dimension

Train_DAT = Train_DAT./( repmat(sqrt(sum(Train_DAT.*Train_DAT)), [size(Train_DAT,1),1]) );
Test_DAT = Test_DAT./( repmat(sqrt(sum(Test_DAT.*Test_DAT)), [size(Test_DAT,1),1]) );

lambda1 = 0.05; % sparse
lambda2 = 0.05; % mean
gamma1  = 10;   % parameter of Eq.(1)
gamma2  = 1;    % parameter of Eq.(1)

% jointly learn the dictionary and the dimension reduction matrix
MaxIter  = 5;
[P,D,C] = JDDLDR(Train_DAT,trls,nDim,numcomps,lambda1,lambda2,gamma1,gamma2,MaxIter);

% prepare the dictionary and testing set
tr_dat  = []; trls = [];
for i = 1:size(D,2)
   tr_dat = [tr_dat D(i).M];
   trls   = [trls repmat(i,[1 size(D(i).M,2)])];
end

tt_dat  =  P'*single(Test_DAT);  clear Test_DAT;
tr_dat  =  tr_dat./( repmat(sqrt(sum(tr_dat.*tr_dat)), [nDim,1]) );
tt_dat  =  tt_dat./( repmat(sqrt(sum(tt_dat.*tt_dat)), [nDim,1]) );

% do collaborative representation based classification
lambda =  0.001;
correct_rate = Fun_CRC(tr_dat,trls,tt_dat,ttls,lambda);

% fid = fopen(['result\demo_JDDLR_result_FRGC.txt'],'a');
% fprintf(fid,'\n%s\n','==========================================');
% fprintf(fid,'%s%8f\n','reco_rate1 = ',correct_rate);
% fclose(fid);