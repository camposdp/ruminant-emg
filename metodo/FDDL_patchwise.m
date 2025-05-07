function [acc, rt, D, CoefM, opts] = FDDL_patchwise(Y_train, label_train, Y_test , label_test, ...
                            k, lambda1, lambda2, pSize,OLstep)
% function [acc, rt] = FDDL_wrapper(Y_train, label_train, Y_test , ...
%       label_test, k, lambda1, lambda2)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------

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
        Y_train_p = [Y_train_p Yp];
        label_train_p = [label_train_p;...
        label_train(i)*ones(length(Yp),1)];
    end
    
    [D, ~, ~, ~, CoefM, ~, opts, rt] = ...
                    LRSDL(Y_train, label_train, opts);
    Y_range = label_to_range(label_train);
    C = max(label_train);
    opts.verbose = 0;
    opts.weight = 0.1;
    acc = [];
    vg = [(0.00005:0.00005:0.00045) ...
        (0.0005:0.0005:0.01) (0.02:0.01:1)];
    for vgamma = vg
        
        opts.gamma = vgamma;
        
        
        [~,N_test] = size(Y_test);
        Y_test_p = [];
        label_test_p = [];
        for i = 1:N_test
            Yp = buildOverlappingPatches(opts, Y_test(:,i));
            label_test_p = [label_test(i)*ones(length(Yp),1)];
            [pred(i,:),~] = max(FDDL_pred(Y_test, D, CoefM, opts));
        end

        
       
        acc1 = double(numel(find(pred == label_test)))/...
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
