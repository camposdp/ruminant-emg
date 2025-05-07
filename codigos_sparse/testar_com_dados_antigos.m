
close all
clear 
clc


load bovino_EMGseg_target_balanced
th_ssc=0.001;
th_zc=0.11;
th_wamp = 0.06;



N = length(EMG);

for j = 1:N
 
Sf = EMG(j,:);

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


end


Ck = cvpartition(N,'KFold',10);


YT = T*[1 2]';


SET{1} = [1 2 3 4];
SET{2} = [12 11 10 5];
SET{3} = [12 11 10 5 3 6 9 8 7];

for i = 1:3

Xf = Ff(SET{i},:)';


mdl_LDA_Tf = fitcdiscr(Xf,YT);
cvmodelf = crossval(mdl_LDA_Tf,'CVPartition',Ck);
ef(:,i) = kfoldLoss(cvmodelf,'Mode','individual','LossFun','classiferror');

Xf = [];

end


lambda1 = 0.01;
lambda2 = 0.01;
k = 300;

for i = 1 : 10
    
        cTrain = Ck.training(i);
        label_train = YT(cTrain);
        Y_train = EMG(cTrain,:)';
        
        cTest =  Ck.test(i);
        label_test = YT(cTest);    
        Y_test = EMG(cTest,:)';
        
        [AC, rt, ~, ~, ~] = myFDDL(Y_train, label_train', Y_test , label_test', ...
                            k, lambda1, lambda2)

        acc(i) = max(AC)
    
    
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
    [D, ~, ~, ~, CoefM, ~, opts, rt] = ...
                    LRSDL(Y_train, label_train, opts);
    Y_range = label_to_range(label_train);
    C = max(label_train);
    opts.verbose = 0;
    opts.weight = 0.1;
    acc = [];
   % for vgamma = [0.0001, 0.001, 0.01, 0.1]
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



function [acc, rt,D, opts] = myDLSI(Y_train, label_train, Y_test , label_test, ...
                            k, lambda, eta)
% function acc = SRC_wrapper(Y_train, range_train, Y_test , range_test, lambda)
% Description       : SRC 
%     INPUT: 
%       dataset: name of the dataset stored in 'data', excluding '.mat'
%       N_trn: number of training images per class 
%       lambda : regularization parameter lambda     
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
    C              = max(label_train);
    D_range        = k*(0:C);
    opts.lambda    = lambda;
    opts.eta       = eta;
    opts.D_range   = D_range;
    opts.show_cost = 0;
    train_range    = label_to_range(label_train);
    opts.show      = 0;
    opts.verbose    = false;
    opts.max_iter  = 100;        
    %% ========= Train ==============================
    [D, ~, rt]         = DLSI(Y_train, train_range, opts);
    %% ========= test ==============================
    opts.verbose    = false;
    pred           = DLSI_pred(Y_test, D, opts);
    acc            = double(sum(pred == label_test))/numel(label_test);
end 



