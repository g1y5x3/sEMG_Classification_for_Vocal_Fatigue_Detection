function [ feature_all, feature_label ] = get_Features_All( features, ...
                                                            DATA, ...
                                                            LABEL, ...
                                                            NOISE, ...
                                                            signatures, ...
                                                            num_channels, ...
                                                            num_signature)

num_gesture = size(DATA,1);

if sum(features(1:7)) > 0
    % Extract time domain features
    feature_TD = cell(num_gesture,num_channels);
    for i = 1:num_gesture
        parfor j = 1:num_channels
            % Raw data
            feature_TD{i,j} = get_TD_features_noise(DATA{i,j},features,NOISE(j));
        end  
    end
end

if features(8) == 1
    % Extract GUSSS ratios
    feature_GR = cell(num_gesture,num_channels);
    parfor i = 1:num_channels
        feature_GR(:,i) = get_GR_features(signatures{1,i}, DATA(:,i), num_signature);
    end
end

% All features are selected, simply concatnate them together
if features(8) == 1 && sum(features(1:7)) > 0

    feature_all = [];
    for j = 1:num_gesture
        feature_tmp = [];
        for i = 1:num_channels
            feature_tmp = horzcat(feature_tmp, feature_TD{j,i});
            feature_tmp = horzcat(feature_tmp, feature_GR{j,i});            
        end
        feature_all = vertcat(feature_all,feature_tmp);    
    end

elseif features(8) == 1 && sum(features(1:7)) == 0

    feature_all = [];
%     for j = 1:num_signature
        feature_tmp = [];
        for i = 1:num_channels
            feature_tmp = horzcat(feature_tmp, feature_GR{1,i});            
        end
%         feature_all = vertcat(feature_all,feature_tmp);    
          feature_all = feature_tmp;
%     end

else

    feature_all = [];
    for j = 1:num_gesture
        feature_tmp = [];
        for i = 1:num_channels
            feature_tmp = horzcat(feature_tmp, feature_TD{j,i});
        end
        feature_all = vertcat(feature_all,feature_tmp);    
    end

end

% Create one column vector of feature labels
feature_label = [];    
if num_gesture > 1
    for j = 1:num_gesture
        feature_label = vertcat(feature_label,LABEL{j});    
    end
else
    feature_label = vertcat(feature_label,LABEL);
end

