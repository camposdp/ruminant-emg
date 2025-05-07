%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daniel Prado de Campos - UTFPR/CPGEI
% 24/05/2019
% DFDL Wrapper
% Single-channel two-classes sEMG DFDL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [acc] = DFDLwrapper(Train1,Train2,Test1,Test2,pars,roc)

%% DFDL
%============ INPUTS ===================%
% pars.numPatches: Number of patches (If patches are chosen randomly)
% pars.overlapStep: Fraction of patch to overlap
%                   e.g. p = 0.5 (50%)
% pars.patchSize: Patch size
% pars.dictsize: Dictionary size per class
% pars.rho: Intra-class regularization term
% pars.max_iter: number of maximum iterations in the main DFDL
% pars.gamma: lambda for ODL
%             l1-norm regularization term
%             for sparse coding in the classification scheme
% pars.lambda: gamma for SRC
%              l1-norm regularization term
%              for sparse coding for estimating L
% pars.paramOMP:
%  OMP pars for mexOMP - see SPAMS document for info
%  http://spams-devel.gforge.inria.fr/doc/html/index.html
%
% Training test:
% Train1: class 1
% Train2: class 2
%
% Training test:
% Test1: class 1
% Test2: class 2
%
% roc: if true, plot ROC curve
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Step 1.1. Building patches

fprintf('Step 1.1. Building patches...\n');

[X, label]  = buildPatches(Train1,Train2,pars);
Y = normc(double(X)); 	% normalize training patches before training
%Y=double(X);
pars.C = [sum(label == 1) sum(label == 2)]; 


%%  Step 1.2 Training dictionaries
fprintf('Step 1.2. Training dictionaries...\n');

[Model, pars] = myDFDL(Y, pars);	% Run the main DFDL
%D1    = Model.Dict(:,:,1);
%D2    = Model.Dict(:,:,2);

%% Step 2.1. Finding proportion between classes
fprintf('Step 2.1 Finding proportion between classes   \n')

[~,Ntrain1] = size(Train1);
[~,Ntrain2] = size(Train2);

ftr1 = zeros(1, Ntrain1);
ftr2 = zeros(1, Ntrain2);
trainlabel = [ones(size(ftr1)) 2*ones(size(ftr2))];

fprintf('Class 1...\n');
for i = 1: Ntrain1

    feature = DFDLonImage(Model, pars, Train1(:,i));
    fprintf('id = %3d, feature = %6f\n', ...
        i, feature);
    ftr1(i) = feature;
end

fprintf('Class 2...\n');
for i = 1: Ntrain2

    feature = DFDLonImage(Model, pars, Train2(:,i));
    fprintf('id = %3d, feature = %6f\n', ...
        i, feature);
    ftr2(i) = feature;
end
    
%%  Step 2.2. Defining threshold (theta)
    fprintf('Step 2.2 Defining threshold (theta) \n')
	F               = [ftr1 ftr2];
	[thresh, signH] = findThrsh(F, trainlabel);
	fprintf('Threshold = %f\n', thresh);    

%%  Step 3. Testing
	fprintf('Step 3. Testing')
    
    [~,Ntest1] = size(Test1);
    [~,Ntest2] = size(Test2);
    
    
	fprintf('Class 1 - test...\n');
	ftest1 = zeros(1, Ntest1);
	ftest2 = zeros(1, Ntest2);
	for i = 1: Ntest1
		
		feature   = DFDLonImage(Model, pars,Test1(:,i));
		ftest1(i) = feature;
		pred1(i)  = -0.5*signH*(2*(feature > thresh) -1) + 1.5;
		fprintf('id = %3d, class = %d\n', ...
			i,  pred1(i));
	end
	fprintf('Class 2 - test...\n');
	for i = 1:  Ntest2
        
		feature   = DFDLonImage(Model, pars, Test2(:,i));
		ftest2(i) = feature;
		pred2(i)  = -0.5*signH*(2*(feature > thresh) -1) + 1.5;
		fprintf('id = %3d, class = %d\n', ...
			i,  pred2(i));
	end
	%% ========= Report results ==============================
	acc1 = sum(pred1 == 1)/numel(pred1);
	acc2 = sum(pred2 == 2)/numel(pred2);
	acc = (sum(pred1 == 1) + sum(pred2 == 2))/numel([pred1 pred2]);
	fprintf('Accuracy: \n');
	fprintf('--------------------------------------\n');
	fprintf('| Class 1    | Class 2   | Overall   |\n')
	fprintf('--------------------------------------\n');
	fprintf('| %4f  | %4f | %4f |\n', 100*acc1, 100*acc2, 100*acc);
	fprintf('--------------------------------------\n');    
    
    
    %% ========= ROC curve ==============================
    if roc == true
	[FAR, MR] = DFDL_ROC(ftest1, ftest2);
	figure;
	plot(FAR, MR, 'bx-'); axis equal;
	hold on;
	plot(1-acc1, 1 - acc2, 'xr');
	title('Receiver Operating Characteristic curve');
	xlabel('Probability of false alarm');
	ylabel('Probability of miss');
	axis([0 1 0 1]);
    end
	% pars

    
end

%% ==============        Functions       =========================%%
    
%%  Build randon patches 
function [X, label] = buildPatches(SigMatrix1,SigMatrix2,pars)
	fprintf('Class 1...');
	[X1] = buildPatches_class(SigMatrix1,pars);
	fprintf('done\nClass 2...');
	[X2] = buildPatches_class(SigMatrix2,pars);
	fprintf('done\n');
	X = [X1 X2]; X = double(X);
	label = [ones(1, size(X1,2)) 2*ones(1, size(X2,2))];
end
%%  Build randon patches for each class
function [X] = buildPatches_class(SigMatrix,pars)  

   
    %Signal matrix N x M:
    %M signals of dimension N
    %M columns, N rows
    S = size(SigMatrix);
    N = S(1);
    M = S(2);
    
	numSignals     = M;
    
	patchSize = pars.patchSize;
    %patchesPerSignal = round(pars.numPatches/numSignals);
    patchesPerSignal = round(N/patchSize);
    
	%Npatches       = patchesPerSignal*numSignals ;
    %X = zeros(patchSize, Npatches);
    %idpatch        = 1;
    X = [];
    for i = 1: numSignals
        data = SigMatrix(:,i);
%         for j = 1: patchesPerSignal
%             
%             xi = (j-1)*patchSize;
%             wi = xi+1: xi+patchSize;
%             
%             
%             data_patch = data(wi);
%             X(:,idpatch) = data_patch(:);
%             idpatch = idpatch + 1;
%         end
    Xi = buildOverlappingPatches(pars, data);
    X = [X Xi];
    end
	
end


%%  Build patches 
function feature = DFDLonImage(Model, pars, signal)
	%  1. Build non-overlapping patches 
	%X = buildNonOverlappingPatches(pars, signal);
    X = buildOverlappingPatches(pars, signal);
	%Y = X ;
    Y = normc(double(X));
	%  SRC 	
	D1                = Model.Dict(:,:,1);
	D2                = Model.Dict(:,:,2);	
	dictsize          = pars.dictsize;
	paramLasso.lambda = pars.gamma;
	paramLasso.eps    = 1e-5;
	D                 = [D1 D2];
	S                 = mexLasso(Y, D, paramLasso);
	S1                = S(1:dictsize,:);
	S2                = S(dictsize+1:end,:);
	R1                = Y - D1*S1;
	R2                = Y - D2*S2;
	e                 = [	sum(R1.^2); sum(R2.^2)];
	[~, pred]         = min(e);
	feature           = sum(pred == 1)/numel(pred);
end
%% Build patches for each signal
function X = buildNonOverlappingPatches(pars, signal)
	p       = pars.patchSize;
	W1 = max(size(signal));
    %signal_overlap=signal(floor(p/2):end-floor(p/2));
    %W2 = max(size(signal_overlap));
	F1       = signal(1:p*floor(W1/p));
    %F2        = signal_overlap(1:p*floor(W2/p));
	X1      = im2col(F1, [p 1], 'distinct');
    %X2      = im2col(F2, [p 1], 'distinct');
    %X = [X1 X2];
    X = X1;
end

function X = buildOverlappingPatches(pars, signal)
	p       = pars.patchSize;
    step = floor(pars.overlapStep*p);
    W = max(size(signal));
    
    xi = 1;
    X=[];
    while xi <= (W-p+1)
    X = [X signal(xi:xi+p-1)];
    xi = xi+step;
    end
    
end

%% Find Threshold
function [thresh, signH] = findThrsh(X_train, trainLabel)
    % find the best thresh for X_train with two classes, 1 dimension
    % return the best thresh and sign of samples from class 1
    % trainLabel = 1 or -1;

    thr_min = min(X_train);
    thr_max = max(X_train);
    step = (thr_max - thr_min)/1000;
    acc = [];
    k = 0;
    h = [];
    for thr = thr_min + step : step:  thr_max - step
        k = k + 1;
        h1 = sum(X_train(find(trainLabel == 1)) > thr) - sum(X_train(find(trainLabel ~=1)) > thr);
        h1 = h1/numel(X_train);
        h = [h h1];
    end
    [accmax, idmax] = max(abs(h));
    maxindex = find(abs(h) == accmax);
    medmax = median(maxindex);

    [~,iid] = min(abs(maxindex - medmax));
    idmax = maxindex(iid);

    thresh = thr_min + (medmax)*step;

    signH = sign(h(idmax));
end


%% construct the ROC 
function [FAR, MR] = DFDL_ROC(ftest1, ftest2)
	N1      = numel(ftest1);
	N2      = numel(ftest2);
	Npoints = 100;
	MR      = zeros(Npoints, 1); % miss rate
	F       = zeros(Npoints, 1); % false alarm rate
	for i = 1:Npoints
		th     = (i)/Npoints;
		acc1   = sum(ftest1 < th)/N1;
		acc2   = sum(ftest2 >= th)/N2;
		MR(i)  = 1 - acc2;
		FAR(i) = 1 - acc1;
	end
end



