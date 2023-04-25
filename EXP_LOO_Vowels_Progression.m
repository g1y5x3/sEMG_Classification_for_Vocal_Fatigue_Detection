%% Inter-Subject Leave-One-Out Classification on Vowels
% Vocally healthy subjects vs. Vocally fatigued subjects
% 1.Assign the class based on their VFI dynamically
% 2.Subject can also be provided as an special case

% function [] = EXP_Inter_H_vs_F(feat_sel, norm)
%% Notes
% 1. This approach assumes there are totally N subjects, the classifier 
%    is trained with N-1 subjects' data and test on the one subject left to
%    compute the diagonsis performance.
clear
clc 

%%
norm = 0;        % Normalization option
ch_sel ='5';     % Channel selection

%% Include library
addpath(genpath('scripts_tools'))
addpath(genpath('EMG Feature Extraction Toolbox'))

%% Initialization Options
% diary EXP_Inter_H_vs_F.out

% Test title
test_title = 'N_vs_F_MATCH';
fprintf('%s\n',test_title)

% Total subjects
%  sub_all = [ 41 25  ...
%              11 19 ];

%  sub_all = [ 41 25 23 34 52 22 04 40 37 15 42 61 44 62 ...
%              11 19 21 43 45 46 47 49 50 53 56 57 59 63];

 sub_all = [ 41 25 23 34 52 22 04 40 37 15 42 61 44 ...
             11 19 21 43 45 46 47 49 50 53 56 57 59 ];

 
% sub_all = 1:61;

% The option of whether duplicate the channel 1 with channel 2
% option_dup = 0;

% Feature selection
% MAV ZC SSC WL WA RMS AR MFL VAR GR
% feat_sel = '000000010';
%feat_sel = '1111111111';
 feat_sel = '1111111001';

feat_txt = {'MAV', 'ZC', 'SSC', 'WL', 'WA', 'RMS', 'AR', 'MFL', 'VAR', 'GR'};
feat_vec = feat_sel - '0';

ch     = ch_sel - '0';

% Partition subjects into healthy group and fatigued group
% sub_h - Vocally Healthy subjects
% sub_f - Vocally Fatigued subjects
[sub_h, sub_f, ~] = partition_subjects(sub_all, 11, 11);

% Initialize parameters 
num_sub = size(sub_all, 2);     % total number of subjects
num_h   = size(sub_h,   2);     % number of vocally healthy subjects
num_f   = size(sub_f,   2);     % number of vocally fatigued subjects

num_cls = 2;            % total number of classes
num_ch  = size(ch,2);   % total number of channels

% Random shuffle the subjects
% subject_healthy = subject_healthy(randperm(length(subject_healthy)));
% subject_dysfunc = subject_dysfunc(randperm(length(subject_dysfunc)));

fprintf('Finish Initialization!\n');
fprintf('Healthy Subjects: %d\n', sub_h)
fprintf('Fatigued Subjects: %d\n', sub_f)
fprintf('Features: %s\n', feat_txt{logical(feat_vec)});
fprintf('Channels: %s\n', ch_sel);

fprintf('Total Subjects: %d\n', num_sub);
fprintf('# of Healthy Subjects: %d\n', num_h);
fprintf('# of Fatigued Subjects: %d\n', num_f);

%% Load data from subjects
fprintf('\nLoading data...\n')

DATA  = cell(num_sub, 5);
LABEL = cell(num_sub, 1);
NOISE = zeros(num_sub, 5);
MVC   = zeros(num_sub, 4);
L     = cell(num_sub, 1);

for i = 1 : num_sub   

    % 1) load subject name
    s = sub_all(i);    
    name = sprintf('R%03d', s);
    disp(name)

    if exist(sprintf('Data/%s',name), 'dir')
        % 2) Assign the data label
        if ismember(s, sub_h)
            label = -1;
        elseif ismember(s, sub_f)
            label = 1;
        end

        % 3) load data samles, labels, signal noise, and mvc (normalization
        % scale)
%         fprintf('2s Noise Pad\n');        
%        [DATA(i,:), LABEL{i}, NOISE(i,:), MVC(i,:)] = loadNorVowels(name, label, 'extracted_voice_pad');
%         fprintf('2s Circular Buffer\n');
%        [DATA(i,:), LABEL{i}, NOISE(i,:), MVC(i,:)] = loadNorVowels(name, label, 'extracted_vowel_buf');
%         fprintf('1s Sliding Window\n');
%         [DATA(i,:), LABEL{i}, NOISE(i,:), MVC(i,:)] = loadVowelSlide1(name, 4000, 1000, false, label);
        fprintf('Fix window\n');
        phrase = {'a_normal', 'u_normal', 'i_normal'};            
        [DATA(i,:), LABEL{i}, NOISE(i,:), MVC(i,:)] = loadVowelFix_audio(name, 4000, label, phrase);    
    end

end

DATA  = DATA(:, ch);
NOISE = NOISE(:, ch);
MVC   = MVC(:, [1,2,3,4]);

fprintf('Finished loading data!\n')

%% Signal Normalization
if norm == 1
    fprintf('Normalize sEMG signals!\n')
    for i = 1 : num_sub
        for c = 1 : num_ch
            % Normalize the original Signal
            DATA{i,c} = DATA{i,c} / MVC(i,c);
            
%             % Normalize the noise
            NOISE(i,c) = NOISE(i,c) / MVC(i,c);
        
        end
    end
end

%% Learning the base signatures for GUSSS average all
% first row  - Vocally Fatigued (Positive)
% second row - Vocally Healthy (Negative)
% tic
% 
% fprintf('Compute GUSSS base signatures ...\n');
signatures = get_Signatures_AVG(DATA, num_cls, num_ch, sub_all, sub_h, sub_f);
% signatures = get_Signatures_ICA(DATA, num_cls, num_ch, sub_all, sub_h, sub_f);
% signatures = get_Signatures_GUSSS(DATA, num_cls, num_ch, sub_all, sub_h, sub_f);
% 
% toc

%% Learning the base signatures for GUSSS average except the left out one
% first row  - Vocally Fatigued (Positive)
% second row - Vocally Healthy (Negative)

% tic
% 
%signatures = cell(num_sub, 1);
%for i = 1 : num_sub
%   fprintf('Compute GUSSS base signatures ...\n');
%   index = 1 : num_sub == i;
%   signatures{i} = get_Signatures_AVG(DATA(~index,:), num_cls, num_ch, ...
%                                      sub_all(:,~index), sub_h, sub_f);
%
%end
% 
% toc


%% Feature extraction
fprintf('Compute features ...\n');
FEAT = cell(num_sub,1);

for i = 1 : num_sub
    s = sub_all(i);
%      [FEAT{i},~] = get_Features_All(feat_vec, DATA(i,:), LABEL{i}, ...
%                                     NOISE(i,:), signatures, num_ch, num_cls);    
%    [FEAT{i},~] = get_Features_All(feat_vec, DATA(i,:), LABEL{i}, ...
%                                   NOISE(i,:), signatures{i}, num_ch, num_cls);
    [FEAT{i},~] = get_Features_New(feat_vec, DATA(i,:), LABEL{i}, ...
                                   NOISE(i,:), signatures, num_ch, num_cls);
%   [FEAT{i},~] = get_Features_New(feat_vec, DATA(i,:), LABEL{i}, ...
%                                  NOISE(i,:), signatures{i}, num_ch, num_cls);
end
fprintf('Finished computing features ...\n');

 %% Post process GUSSS ratio
% 
% % FEAT = FEAT_COPY;
% % 
% for i = 1 : num_sub
% %     FEAT{i} = FEAT{i}(:,[2,4,6,8]);
%     FEAT{i} = trapmf(FEAT{i},[0 300 1000000 1000000]);
% end

%% Normalize the feature vectors
[Z, mu, sigma] = zscore(cell2mat(FEAT));

FEAT_N = cell(num_sub, 1);
for s = 1 : num_sub
    for d =  1 : size(mu,2)
        FEAT_N{s,1}(:,d) = (FEAT{s}(:,d) - mu(d)) / sigma(d);        
    end    
end

FEAT = FEAT_N;

%% Leave-One-Out Crossvalidation
for lo_pen = 1 : 0.1 : 1
size_train          = zeros(num_sub, num_cls);
size_test           = zeros(num_sub, num_cls);
size_total          = zeros(num_sub, num_cls);
accuracy_train      = zeros(num_sub, 1);
accuracy_test       = zeros(num_sub, 1);
accuracy_validation = zeros(num_sub, 1);

fprintf('\nLeft Out percentage: %.2f\n', lo_pen);

% For VFI correlation analysis
positive_per        = zeros(num_sub,1);
VFI1_test           = zeros(num_sub,1);
load('Notes/VFI.mat')

for i = 1 : num_sub    
    if exist(sprintf('Data/R%03d', sub_all(i)), 'dir')

    % Partition the testing data
    fprintf('\nTesting Subject: R%03d\n', sub_all(i));
    VFI1_test(i) = VFI1(sub_all(i));
    fprintf('VFI-1 Score: %d\n', VFI1(sub_all(i)));    

    if lo_pen ~= 0
        % Assign the testing data (from the current subject)
        feature_lo = FEAT{i};    
        label_lo   = LABEL{i};

        size_test(i,1) = sum(label_lo == 1);
        size_test(i,2) = sum(label_lo == -1);

        N = size(label_lo, 1);
        num_lo = round(N * lo_pen);
        x = randsample(N,num_lo);

        feature_test = feature_lo(x,:);
        label_test   = label_lo(x,:);

        y = (1 : N)';
        y(x) = [];

        % Assign the training data (from the left subjects)
        FEAT_TMP  = FEAT;
        LABEL_TMP = LABEL;    
        FEAT_TMP(i,:)  = [];
        LABEL_TMP(i,:) = [];

        FEAT_TMP = FEAT_TMP(~cellfun('isempty',FEAT_TMP));
        LABEL_TMP = LABEL_TMP(~cellfun('isempty',LABEL_TMP));    

        feature_train = cell2mat(FEAT_TMP);
        label_train   = cell2mat(LABEL_TMP);
        feature_train = vertcat(feature_train, feature_lo(y,:));
        label_train = vertcat(label_train, label_lo(y,:));

    else
        feature_train = cell2mat(FEAT);
        label_train   = cell2mat(LABEL);
        feature_test = FEAT{i};    
        label_test   = LABEL{i};        
        
    end    
%     display(size(feature_test,1))

    size_train(i,1) = sum(label_train == 1);
    size_train(i,2) = sum(label_train == -1); 
   
    % The total amount of training + testing data
    size_total(i,:) = size_test(i,:) + size_train(i,:);
 
    % Calculate the ratio between the positive and negative classes
    ratio = size_train(i,2) / size_train(i,1);
    fprintf('Ratio: %f\n', ratio);
   
    % Train the SVM classifier and perform cross-validation
    fprintf('Training...\n');
     SVMModel = fitcsvm(feature_train,label_train,...
                       'Standardize',true,...                          
                       'KernelFunction','RBF',...
                       'ClassNames',[1,-1],...
                       'KernelScale', 1,...
                       'Cost',[0,ratio;1,0],...
                       'OutlierFraction',0.05);
%    % Some of the parameters that have produced good results
%                      'Standardize',true,...                  
%                       'KernelScale', 0.7e01,...                      
%                      'KernelScale', 3.25e01,...
%                      'BoxConstraint', 1e-1,...

%   c = cvpartition(sum(size_train(i,:)),'KFold',10);
% 
%   opts = struct('Optimizer','bayesopt',...
%                 'ShowPlots',false,...
%                 'Verbose', 0,...
%                 'CVPartition',c,...
%                 'AcquisitionFunctionName','expected-improvement-plus');
%   
%   svmmod = fitcsvm(feature_train,label_train,...
%                    'Cost',[0,ratio;1,0],...
%                    'KernelFunction','rbf',...
%                    'OptimizeHyperparameters','auto',...
%                    'HyperparameterOptimizationOptions',opts);    
% 
%                    %'Standardize',true,...        
% 
%   SVMModel = fitcsvm(feature_train,label_train,...
%                      'Cost',[0,ratio;1,0],...
%                      'KernelFunction','rbf',...
%                      'BoxConstraint',svmmod.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint,...
%                      'KernelScale',svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelScale);    

    % Compute the error of training
    accuracy_train(i) = 1 - resubLoss(SVMModel);
    fprintf('Training Accuracy: %f\n', accuracy_train(i));
                 
    % Cross validate the SVM classifier and compute the validation error
    CVSVMModel = crossval(SVMModel);
    accuracy_validation(i) = 1 - kfoldLoss(CVSVMModel);
    fprintf('Validation Accuracy: %f\n', accuracy_validation(i));
    
    % Apply the classifier on the testing dataset and compute final error
    positive_per_fold = zeros(10,1);
    
    % 10 - fold cross-validation
    for f = 1 : 10
        CompactSVMModel = CVSVMModel.Trained{f};
        loss_tmp = loss(CompactSVMModel, feature_test, label_test);
        accuracy_test(i) = accuracy_test(i) + (1 - loss_tmp);
        
        % Measure the amount of predictive positive samples
        label_predict = predict(CompactSVMModel, feature_test);        
        positive_per_fold(f) = sum(label_predict == 1)/size(label_predict,1);                  
    end
    accuracy_test(i) = accuracy_test(i)/10;
    positive_per(i)   = mean(positive_per_fold);
    
    fprintf('Testing Accuracy: %f\n', accuracy_test(i));
    fprintf('Percentage of Positive Samples: %f\n', positive_per(i));    
    end
end

accuracy_dysfunc = accuracy_test(logical(sum(sub_all==sub_f',1)));
accuracy_healthy = accuracy_test(logical(sum(sub_all==sub_h',1)));

% Display the numerical results
fprintf('\nAverage Training Accuracy: %.2f%%\n', mean(accuracy_train(accuracy_train ~= 0))*100);
fprintf('Average Validation Accuracy: %.2f%%\n', mean(mean(accuracy_validation(accuracy_train ~= 0)))*100);
fprintf('Average Testing Accuracy: %.2f%%\n', mean([accuracy_dysfunc;accuracy_healthy])*100);
fprintf('Sensitivity(True Positive): %.4f\n',mean(accuracy_dysfunc));
fprintf('Specificity(True Negative): %.4f\n',mean(accuracy_healthy));
%disp('Training')
%disp(accuracy_train)
%disp('Validation')
%disp(accuracy_validation)
%disp('Testing')
%disp(accuracy_test)

end
% diary off

%% Save the results
% Save the classification results
% save(sprintf('InterResults_Match/%s.mat', test_title),...
%      'feature_vector', ...
%      'accuracy_train', ...
%      'accuracy_validation', ...
%      'accuracy_test', ...
%      'accuracy_dysfunc', ...
%      'accuracy_healthy', ...
%      'size_train',...
%      'size_test', ...
%      'size_total', ...
%      'subject_all', ...
%      'num_subjects',...
%      'subject_healthy', ...
%      'num_healthy', ...
%      'subject_dysfunc', ...
%      'num_dysfunc',...
%      'positive_per', ...
%      'VFI1_test');
% 
% fprintf('Finish saving the results and individual classifier!\n');

% %% Correlation Analysis
% index = accuracy_test ~= 0;
% positive_per_tmp = positive_per(index);
% VFI1_test_tmp = VFI1_test(index);
% 
% rho = corr(positive_per_tmp, VFI1_test_tmp);
% fprintf('Correlation between VFI-1 and Positive Samples Percentage: \n');
% fprintf('%f\n', rho);
% 
% % Calculate the average positive samples percentage for each VFI-1
% positive_per_vfi = zeros(max(VFI1_test_tmp)+1,1);
% for vfi1 = min(VFI1_test_tmp) : max(VFI1_test_tmp)
%     positive_per_vfi(vfi1+1) = mean(positive_per_tmp(VFI1_test_tmp == vfi1));
% end
% 
% bar((0:28),positive_per_vfi)
% xlim([-1 30])
% ylim([0 1])
% title(sprintf('Corr = %.3f',rho));
% ylabel('Positive Rate')
% xlabel('VFI-1 for testing subjects')
% % saveas(gcf,'PR_Corr.png')
% %%
% figure;
% bar((0:28),positive_per_vfi)
% hold on
% bar((0:10),positive_per_vfi(1:11),'c')
% bar((11:28),positive_per_vfi(12:29),'b')
% xlim([-1 30])
% ylim([0 1])
% ylabel('Positive Rate','FontSize',15)
% xlabel('VFI-1 for testing subjects','FontSize',15)
% legend('Vocally Fatigued', 'Vocally Healthy','Location','northwest')
