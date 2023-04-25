function [ feature_vector ] = get_TD_features_noise( signal, features, noise_thres )
% Create a feature vector for the EMG signals, which contains all time
% domain features.
% feature_vector = [MAV ZC SSC WL WA RMS AR]
% MAV   -   Mean absolute value
% ZC    -   Zero crossing
% SSC   -   Slope sign changes
% WL    -   Waveform length
% WA    -   Willison amplitude
% RMS   -   Root Mean Square
% AR    -   Autoregressive coefficients
% --------------------------------------------------
N = size(signal,1);         % total number of the repetitions
L = size(signal,2);         % total length of the signal

%***TO-DO*** HOW TO FIND THE OPTIMAL THRESHOLD VALUE
% noise_thres = 0.001;       % threshold for the noise when counting zc 
                            % and ssc unit in volts
% Mean Absolute Value
if features(1) == 1
    mav = sum(abs(signal),2) / L;
end

% Zero Crossings
if features(2) == 1
    zc = zeros(N,1);
    for rep = 1:N
        for i = 1:L-1
            if ((signal(rep,i) > 0 && signal(rep,i+1) < 0) || ...
                    (signal(rep,i) < 0 && signal(rep,i+1) > 0)) && ...
                    abs(signal(rep,i) - signal(rep,i+1)) >= noise_thres
                % When detect the signal pass across zero, increment zc
                zc(rep) = zc(rep) + 1;
            end
        end
    end
end

% Slope Sign Changes
if features(3) ==1
    ssc = zeros(N,1);
    for rep = 1:N
        for i = 2:L-1
            if ((signal(rep,i) > signal(rep,i-1) && ...
                    signal(rep,i) < signal(rep,i+1)) || ...
                    (signal(rep,i) < signal(rep,i-1) && ...
                    signal(rep,i) > signal(rep,i+1))) && ...
                    (abs(signal(rep,i) - signal(rep,i+1)) >= noise_thres || ...
                    abs(signal(rep,i) - signal(rep,i-1)) >= noise_thres)
                ssc(rep) = ssc(rep) + 1;
            end
        end
    end
end

% Waveform Length
if features(4) == 1
    wl = sum(abs(signal(:,2:L) - signal(:,1:L-1)),2);
end

% Willison Amplitude
if features(5) == 1
    
    wa = sum(abs(signal(:,1:L-1) - signal(:,2:L) > noise_thres),2);

end

% Root Mean Square
if features(6) == 1
    rms = sqrt(sum(signal.^2,2) / L);
end

% Atuoregressive Coefficients
if features(7) == 1
    arcof = zeros(N,4);
    for i = 1 : N
       m = ar(signal(i,:),4);
       arcof(i,:) = m.a(2:5);
    end
end

% Combine all feature vectors
feature_vector = [];
if features(1) == 1
    feature_vector = horzcat(feature_vector,mav);
end

if features(2) == 1
    feature_vector = horzcat(feature_vector,zc);
end

if features(3) == 1
    feature_vector = horzcat(feature_vector,ssc);
end

if features(4) == 1
    feature_vector = horzcat(feature_vector,wl);
end

if features(5) == 1
    feature_vector = horzcat(feature_vector,wa);
end

if features(6) == 1
    feature_vector = horzcat(feature_vector,rms);
end

if features(7) == 1
    feature_vector = horzcat(feature_vector,arcof);
end

end