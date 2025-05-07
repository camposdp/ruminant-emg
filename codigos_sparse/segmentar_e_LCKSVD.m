%RESULTADOS
% [0.748 0.786 0.815 0.815 0.818 0.828 0.823]


close all
clear 
clc

load EMG_raw_labels
load EMGre_exp

fsa = 2000; %frequência de amostragem





valpha = 0.001;
vbeta = 0.001;

BL = [4.45 4.7]*1e5;

%Wcrit = [0.015*fsa 0.6*fsa];


W = [0.2 0.4]*fsa;  
%W = [.2 .4]*fsa;
Wcrit = [0.10 1]*fsa;
Nb = W(1);
Nf = W(2);

%% 3) Segmentação DTOS e extração de características

%Wsmooth = 0.1*fsa;
%k=3.5;

%k = 200;
%sparsitythres = 5; % sparsity prior


Wsmooth = 0.15*fsa;
Kstd = 4.0;




%         
% sparsitythres = 5;
% d = 1;
% k = 300;
% %k = 300;
% %sparsitythres = 5; % sparsity prior
% 
% lambda = 0.01;
% eta = 0.1;

x = smoothdata(abs(EMGn'),'movmean',Wsmooth);
th = mean(EMGt_re(BL(1):BL(2)))+Kstd*std(EMGt_re(BL(1):BL(2)));


% Usar pra obter os segmentos para extração posterior 
[xi,xf,wT] = DTSeg(EMGn',Wcrit,x,th,T);

xi = xi(1:end-1);
xf = xf(1:end-1);
wT = wT(1:end-1);



Clab{1} = find(wT==1);
Clab{2} = find(wT==2);
Clab{3} = find(wT==3);
K = 0.8;

    wTv = [];
    xiv = [];
    wTt = [];
    xit = [];

for i = 1:2
    %N = length(Clab{i});
    N = 1500;
    Nval = randperm(N,round(N*K));
    Ntest = 1:N; 
    Ntest(Nval) = [];
    
    wTv = [wTv wT(Clab{i}(Nval))];
    xiv = [xiv xi(Clab{i}(Nval))];
    wTt = [wTt wT(Clab{i}(Ntest))];
    xit = [xit xi(Clab{i}(Ntest))];
end


label_test = wTt;
label_train = wTv;


N_t = length(label_test);
N_v = length(label_train);

%% 1) Preparar Matriz de EMG

w0 = W(1)+W(2);
d=1;
w = w0/d;
Y_train = zeros(w,N_v);
Y_test = zeros(w,N_t);

for i = 1:N_v
Y_train(:,i) = decimate(EMGn(xiv(i)-Nb:xiv(i)+Nf-1),d); 
end

for i = 1:N_t
Y_test(:,i) = decimate(EMGn(xit(i)-Nb:xit(i)+Nf-1),d); 
end




% %% 3) Parâmetros LC-kSVD
% 
% 
% 
% 
%  acc = zeros(1, 2);
% rt = zeros(1, 2);
%     C = max(label_train);
%     H_train = lcksvd_buildH(label_train);
%     H_test = lcksvd_buildH(label_test);
% 
% 
%     sqrt_alpha = sqrt(valpha); % weights for label constraint term
%     sqrt_beta = sqrt(vbeta); % weights for classification err term
%     dictsize = C*k; % dictionary size
%     iterations = 50; % iteration number
%     iterations4ini = 20; % iteration number for initialization
% 
% %% 4) Aprendizado do dicionário 
%     % get initial dictionary Dinit and Winit
%     fprintf('\nLC-KSVD initialization... ');
%     [Dinit,Tinit,Winit,Q_train] = initialization4LCKSVD(Y_train, H_train, ...
%         dictsize, iterations4ini, sparsitythres);
%     fprintf('done!');
%     
% %% 5) ========= LCKSVD1 ==============================  
%     % run LC K-SVD Training (reconstruction err + class penalty)
%     fprintf('\nDictionary and classifier learning by LC-KSVD1...');
%     tic;
%     [D1,X1,T1,W1] = labelconsistentksvd1(Y_train, Dinit, Q_train, Tinit, ...
%         H_train,iterations,sparsitythres,sqrt_alpha);
%     rt(1) = toc;
%     fprintf('done!');
%     % classification process
%     [prediction1, acc(1)] = classification(D1, W1, Y_test, H_test, sparsitythres);
%     fprintf('\nFinal recognition rate for LC-KSVD1 is : %.03f \n', acc(1));
% %% 6) ========= LCKSVD2 ==============================  
%     % run LC k-svd training (reconstruction err + class penalty + classifier err)
%     fprintf('\nDictionary and classifier learning by LC-KSVD2...');
%     tic;
%     [D2,X2,T2,W2] = labelconsistentksvd2(Y_train, Dinit, Q_train, ...
%         Tinit, H_train, Winit, iterations, sparsitythres, sqrt_alpha, sqrt_beta);
%     rt(2) = toc;
%     fprintf('done!\n');
%     %classification
%     [prediction2, acc(2)] = classification(D2, W2, Y_test, H_test, sparsitythres);
%     fprintf('\nFinal recognition rate for LC-KSVD2 is : %.03f \n', acc(2));
%         
% %     
% % 
% k = 300;
%  
% lambda = 0.01;
% eta = 1;
% [acc(3), rt] = DLSI_wrapper(Y_train, label_train, Y_test , label_test, ...
%                             k, lambda, eta)


lambda1 = 0.01;
lambda2 = 0.01;

Ki = [100 200 300 400 500 600 700];
for i = 1:length(Ki)
k= Ki(i);

[AC, rt] = FDDL_wrapper(Y_train, label_train, Y_test , label_test, ...
                            k, lambda1, lambda2)

acc(i) = max(AC);

end
                        
                        
                        
% lambda = 0.01;
% 
% range_train = label_to_range(label_train);
% range_test = label_to_range(label_test);

% function acc = SRC_wrapper(Y_train, range_train, Y_test , range_test, lambda)
% Description       : SRC 
%     INPUT: 
%       dataset: name of the dataset stored in 'data', excluding '.mat'
%       N_trn: number of training images per class 
%       lambda : regularization parameter lambda    

% 
% acc(5) = SRC_wrapper(Y_train, range_train, Y_test , range_test, lambda)
%     

% 
% 
% 
% lambda1 = 0.01;
% lambda2 = 0.01;
% alpha = 0.01;
% k = 100;
%     
% [acc_temp, rt] = D2L2R2_wrapper(Y_train, label_train, Y_test, label_test,...
%                         k, lambda1, lambda2, alpha)
% 
% acc(6) = max(acc_temp);                    
% 
% k = 100;
% k0 = 5;
% lambda = 0.001;
% eta = 0.01;    
% 
% 
% 
% [acc(7), rt] = COPAR_wrapper(Y_train, label_train, Y_test , label_test, ...
%                             k, k0, lambda, eta)
% 
% %                         
% lambda1 = 0.001;
% lambda2 = 0.01;
% lambda3 = 0.01;
% k = 300;
% k0 = 5;
% 
% [acc(8), rt] = LRSDL_wrapper(Y_train, label_train, Y_test , label_test, ...
%                             k, k0, lambda1, lambda2, lambda3)
% 
%   bar(acc)                      
%   str = {'LC-kSVD1','LC-kSVD2','DLSI','FDDL','SRC', 'D2L2R2' 'COPAR','LRSDL'}                      
%   set(gca, 'XTickLabel',str, 'XTick',1:numel(str))                      
%   xtickangle(45)                      





 
% 
% 
% [acc, rt] = DLSI_wrapper(Y_train, label_train, Y_test , label_test, ...
%                             k, lambda, eta)

% 
% lambda1 = 0.01;
% lambda2 = 0.01;


% [acc(2,i), rt] = FDDL_wrapper(Y_train, label_train, Y_test , label_test, ...
%                             k, lambda1, lambda2)
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

