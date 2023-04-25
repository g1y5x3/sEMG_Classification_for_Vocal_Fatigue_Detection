function [ signatures ] = get_Signatures_ICA( DATA, num_cls, num_ch, sub_all, ...
                                          sub_h, sub_f)

num_sub = size(sub_all,2);                                      
                                      
% create a cell that contains all selected samples for calculating the base
% signatures
SIG_BASE = cell(num_cls, num_ch);

% Approcah 1: use the average across all subjects
for i = 1 : num_sub
    s = sub_all(i);
    if ismember(s, sub_f)
        sig_ind = 1;
    elseif ismember(s, sub_h)
        sig_ind = 2;
    end
    

    for c = 1 : num_ch
      
        sig_sel = DATA{i,c};
        SIG_BASE{sig_ind,c} = vertcat(SIG_BASE{sig_ind,c}, sig_sel);
%        disp(size(SIG_BASE{sig_ind,c}))
    end
end

signatures = cell(1,num_ch);
for i = 1:num_ch
    for s = 1 : num_cls
        [~, AinitS, ~] = fastica(SIG_BASE{s,i}, 'numOfIC',1);
%         AinitS = fliplr(eye(size(SIG_BASE{s,i},1)));
        [icasig,mixA,~] = fastica(SIG_BASE{s,i},...
                                  'numOfIC',1,...
                                  'approach','defl',...
                                  'verbose','off',...
                                  'g','tanh',...
                                  'initGuess',AinitS,...
                                  'epsilon',0.00001,...
                                  'maxNumIterations',1000);      
        mixAmean = mean(mixA,1);
        signatures{1,i}(s,:) = mixAmean * icasig;
    end
end

% signatures = cell(1,num_channel);
% parfor i = 1:num_channel
%     for s = 1 : num_class
%         signatures{1,i}(s,:) = mean(DATA_TRAIN{s,i},1)/max(abs(mean(DATA_TRAIN{s,i},1)));
%     end
% end

end

