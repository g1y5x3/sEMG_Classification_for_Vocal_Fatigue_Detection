%% Extract EMG samples from the original data based on the voice channel

% *** Using zero crossings threshold method to improve the signal detection
% accuracy

clear,clc, close all

%% Initialization of filepath and subject & gesture information
filepath = '/home/iris/yg5d6/NIH/';
savepath = 'Data';

subject_name = 'R067S1';

gesture_name = {'a_normal',...
                'u_normal',...
                'i_normal'};
            
% gesture_name = {'sentence1','sentence2'};

gesture_num = size(gesture_name, 2);


% Parameters for examine particualr repetition
ges_start = 1;      % between 1 to 6
rep_start = 1;      % between 1 to 55

% Based on Frequency or Amplitude
% 1 - amplitude
% 2 - frequency
mov_ave_base = 1;

% FOR DEBUGGING
% 0: plot and print all
% 1: plot and print only necessary
% 2: plot and print nothing
% plot_count is basically a mod divide to decide how much sample you want
% to inspect before moving on
                                                                                                                    
plot_ctrl = 1;
plot_count = 1;

% % Based on what ever subjects that was mentioned in the excel note sheets
% % from data collection
% discard_samples = {[];...           % a normal
%                    [];...           % u normal
%                    []};             % i normal
%                
% save(sprintf('Data/%s/discard_vowels.mat',subject_name), ...
%      'discard_samples');               

%% Get the noise average from sample
% Start the procedure
disp(subject_name)

% Load the noise data file
% *************************************************************************
% switch to 'baseline1' for subjects after R061
% *************************************************************************
if exist(sprintf('%s/%s/%s.mat',filepath,subject_name,'baseline1'),'file') ~= 0

    load(sprintf('%s/%s/%s.mat',filepath,subject_name,'baseline1'));

    % Get the sampling rate for the voice channel (either 4k or 20k)
    rate_voice = samplerate(5);
    rate_emg = samplerate(1);    % sampling rate for the emg channel

    % The noise for each channel will be based on the average of data collected
    noise_ave = zeros(4,1);
    noise_ave(1) = mean(data(datastart(1):dataend(1)));
    noise_ave(2) = mean(data(datastart(2):dataend(2)));
    noise_ave(3) = mean(data(datastart(3):dataend(3)));
    noise_ave(4) = mean(data(datastart(4):dataend(4)));
    noise_ave_voice  = mean(data(datastart(5):dataend(5)));

    % The total length of each channel is 2s, this is a fixed value
    noise = zeros(4,rate_emg * 2);
    noise(1,:) = data(datastart(1):dataend(1));
    noise(2,:) = data(datastart(2):dataend(2));
    noise(3,:) = data(datastart(3):dataend(3));
    noise(4,:) = data(datastart(4):dataend(4));
    noise_voice = data(datastart(5):dataend(5));

    clear blocktimes com comtext data dataend datastart firstsampleoffset ...
        rangemax rangemin samplerate tickrate titles unittext unittextmap

    %========== PRINT DEBUG ==========
    if plot_ctrl < 2
        fprintf('Completed loading data and noise!\n')
    end
 
else
    if plot_ctrl < 2
        fprintf('Noise file does not exist!\n')
    end
end

%% Load previously stored noise signature
load(sprintf('Data/%s/extracted_voice_%s/statistics.mat',subject_name, 'pad'));

%% Extract EMG activity based on the activation
extract_win = 8000;    

for z = ges_start : gesture_num
    
    if exist(sprintf('%s/%s/%s.mat',filepath,subject_name,gesture_name{z}),'file') ~= 0 
        load(sprintf('%s/%s/%s.mat',filepath,subject_name,gesture_name{z}));
        
        clear blocktimes com comtext firstsampleoffset ...
            rangemax rangemin samplerate tickrate titles unittext unittextmap   

        data_extract = cell(1, 4);
        for i = 1 : 4
            data_extract{i} = zeros(55, extract_win);
        end

        % the repetation of the gesture
%         for rep = rep_start :  55
        rep = rep_start;
        while rep <= 55 
            if rep == 0
                rep = rep_start;
            end

            % THIS IS THE RAW DATA WITHOUT AN   Y MODIFICATON
            % assign each channel the raw data for every repetition
            ch_1     = data(datastart(1,rep):dataend(1,rep)); 
            ch_2     = data(datastart(2,rep):dataend(2,rep)); 
            ch_3     = data(datastart(3,rep):dataend(3,rep)); 
            ch_4     = data(datastart(4,rep):dataend(4,rep)); 
            ch_voice = data(datastart(5,rep):dataend(5,rep));

            index_start      = round(signal_start{z}(rep) * rate_emg);
            index_end        = round(signal_end{z}(rep) * rate_emg);
            extracted_length = round(signal_length{z}(rep) * rate_emg);

%             % Extract the signal
%             extracted_length = index_end - index_start + 1;
% 
%             signal_start{z}(rep) = index_start/rate_emg;
%             signal_end{z}(rep) = index_end/rate_emg;
%             signal_length{z}(rep) = extracted_length/rate_emg;   
            if plot_ctrl < 1 && mod(rep,plot_count) == 0
                %% Plot the signal extraction marker
                figure('units','normalized','outerposition',[0 0 1 1])    
                subplot(5,1,1)
                plot(ch_1)
                ylim([-0.05 0.05])
                hold on
                line([index_start index_start], [-0.05 0.05],'Color','red','LineWidth',3);
                line([index_end index_end], [-0.05 0.05],'Color','black','LineWidth',3);
                hold off            
                title(sprintf('%s Repetition: %d',gesture_name{z}, rep),'Interpreter', 'none')
                subplot(5,1,2)
                plot(ch_2)
                ylim([-0.05 0.05])                    
                hold on
                line([index_start index_start], [-0.05 0.05],'Color','red','LineWidth',3);
                line([index_end index_end], [-0.05 0.05],'Color','black','LineWidth',3);
                hold off            
                subplot(5,1,3)
                plot(ch_3)
                ylim([-0.05 0.05])                    
                hold on
                line([index_start index_start], [-0.05 0.05],'Color','red','LineWidth',3);
                line([index_end index_end], [-0.05 0.05],'Color','black','LineWidth',3);
                hold off            
                subplot(5,1,4)
                plot(ch_4)
                ylim([-0.05 0.05])                       
                hold on
                line([index_start index_start], [-0.05 0.05],'Color','red','LineWidth',3);
                line([index_end index_end], [-0.05 0.05],'Color','black','LineWidth',3);
                hold off            
                subplot(5,1,5)
                plot(ch_voice)
                ylim([-0.1 0.1])

                rep_button = uicontrol('Style', 'pushbutton', ...
                                       'String', 'Repeat', ...
                                       'Position', [745 900 100 25],...
                                       'Callback', 'continue_flag = 0; uiresume(gcbf)');

                cont_button = uicontrol('Style', 'pushbutton', ...
                                        'String', 'Continue', ...
                                        'Position', [945 900 100 25],...
                                        'Callback', ...
                                        'continue_flag = 1; uiresume(gcbf);');

                inc_button = uicontrol('Style', 'pushbutton', ...
                                        'String', 'Increase Threshold', ...
                                        'Position', [300 900 175 25],...
                                        'Callback', 'lower_fact = lower_fact - 0.05;'); 

                dec_button = uicontrol('Style', 'pushbutton', ...
                                        'String', 'Decrease Threshold', ...
                                        'Position', [500 900 175 25],...
                                        'Callback', 'lower_fact = lower_fact + 0.05;');

                amp_button = uicontrol('Style', 'pushbutton', ...
                                        'String', 'Amplitude', ...
                                        'Position', [300 875 175 25],...
                                        'Callback', 'mov_ave_base = 1;');

                fre_button = uicontrol('Style', 'pushbutton', ...
                                        'String', 'Frequency', ...
                                        'Position', [500 875 175 25],...
                                        'Callback', 'mov_ave_base = 2;'); 

                back_button = uicontrol('Style', 'pushbutton', ...
                                        'String', 'Back', ...
                                        'Position', [1145 900 100 25],...
                                        'Callback', ...
                                        'continue_flag = 1; uiresume(gcbf); rep = rep-2;');                                        

                uiwait(gcf);          
            end
            
            % 1. Use the circular buffer approach to create a 8000 length
            for s = 1 : extract_win
%                 display(mod(index_start+s-1,extracted_length)+index_start)
                data_extract_buf{1}(rep,s) = ch_1(mod(s,extracted_length)+index_start);
                data_extract_buf{2}(rep,s) = ch_2(mod(s,extracted_length)+index_start);
                data_extract_buf{3}(rep,s) = ch_3(mod(s,extracted_length)+index_start);
                data_extract_buf{4}(rep,s) = ch_4(mod(s,extracted_length)+index_start);                
            end

            % 2. Use only between the beginning and ending of the signal 
            %    and fill the rest of the extraction window with fixed
            %    noise.
            data_extract_pad{1}(rep,1 : extracted_length) = ch_1(index_start : index_end);
            data_extract_pad{2}(rep,1 : extracted_length) = ch_2(index_start : index_end);
            data_extract_pad{3}(rep,1 : extracted_length) = ch_3(index_start : index_end);
            data_extract_pad{4}(rep,1 : extracted_length) = ch_4(index_start : index_end);

            data_extract_pad{1}(rep,extracted_length+1 : extract_win) = noise(1,extracted_length+1 : extract_win);
            data_extract_pad{2}(rep,extracted_length+1 : extract_win) = noise(2,extracted_length+1 : extract_win);
            data_extract_pad{3}(rep,extracted_length+1 : extract_win) = noise(3,extracted_length+1 : extract_win);
            data_extract_pad{4}(rep,extracted_length+1 : extract_win) = noise(4,extracted_length+1 : extract_win);

            % 3. Use the beginning, shift it forward and fill the rest of
            %    the window with random sampled noise buffer.
            slide_length = 8000 - index_start + 1;
            data_extract_slide{1}(rep,1 : slide_length) = ch_1(index_start : extract_win);
            data_extract_slide{2}(rep,1 : slide_length) = ch_2(index_start : extract_win);
            data_extract_slide{3}(rep,1 : slide_length) = ch_3(index_start : extract_win);
            data_extract_slide{4}(rep,1 : slide_length) = ch_4(index_start : extract_win);
            
            remain_length = 8000 - slide_length;
            noise_index   = randi(extract_win, remain_length, 1);
            data_extract_slide{1}(rep,slide_length+1 : extract_win) = noise(1,noise_index);
            data_extract_slide{2}(rep,slide_length+1 : extract_win) = noise(2,noise_index);
            data_extract_slide{3}(rep,slide_length+1 : extract_win) = noise(3,noise_index);
            data_extract_slide{4}(rep,slide_length+1 : extract_win) = noise(4,noise_index);            

            rep = rep+1;
            
            close all
            
        end

        if plot_ctrl < 2
            fprintf('Completed extracting gestures %s!\n', gesture_name{z})
        end

        % 1. Save the data that uses circular buffer approach
        data   = data_extract_buf;
        length = signal_length{z};

        if exist(sprintf('Data/%s/extracted_vowel_%s', subject_name, 'buf'),'dir') ~= 7
            mkdir(sprintf('Data/%s/extracted_vowel_%s', subject_name, 'buf'))
            
            save(sprintf('Data/%s/extracted_vowel_%s/%s.mat',subject_name, 'buf', gesture_name{z}), ...
                 'data', ...
                 'length');            
        else
            save(sprintf('Data/%s/extracted_vowel_%s/%s.mat',subject_name, 'buf', gesture_name{z}), ...
                 'data', ...
                 'length');
        end

        % 2. Save the data that uses padding approach
        data = data_extract_pad;
        length = signal_length{z};

        if exist(sprintf('Data/%s/extracted_voice_%s', subject_name, 'pad'),'dir') ~= 7
            mkdir(sprintf('Data/%s/extracted_voice_%s', subject_name, 'pad'))
            
            save(sprintf('Data/%s/extracted_voice_%s/%s.mat',subject_name, 'pad', gesture_name{z}), ...
                 'data', ...
                 'length');           
        else
            save(sprintf('Data/%s/extracted_voice_%s/%s.mat',subject_name, 'pad', gesture_name{z}), ...
                 'data', ...
                 'length');
        end
         
        % 3. Save the data that uses sliding and randomly sampled noise
        %    approach
        data = data_extract_slide;
        length = signal_length{z};

        if exist(sprintf('Data/%s/extracted_vowel_%s', subject_name, 'slide'),'dir') ~= 7
            mkdir(sprintf('Data/%s/extracted_vowel_%s', subject_name, 'slide'))
            
            save(sprintf('Data/%s/extracted_vowel_%s/%s.mat',subject_name, 'buf', gesture_name{z}), ...
                 'data', ...
                 'length');            
        else
            save(sprintf('Data/%s/extracted_vowel_%s/%s.mat',subject_name, 'slide', gesture_name{z}), ...
                 'data', ...
                 'length');
        end
        
        if plot_ctrl < 2
            fprintf('Completed saving extracted gestures!\n')
        end
    else
        if plot_ctrl < 2
            fprintf(sprintf('%s file does not exist!\n', gesture_name{z}));
        end
    end
   
end

%% Save the statistics of the sample selection
% Only saves the statics data whenever starting collecting from the first
% gesture and includes all of the repetitions.

% if ges_start == 1 && rep_start == 1
%     if exist('extract_win', 'var') ~= 0
%         save(sprintf('Data/%s/extracted_vowel_%s/statistics.mat',subject_name, 'buf'),...
%              'signal_start','signal_end','signal_length');
%         save(sprintf('Data/%s/extracted_voice_%s/statistics.mat',subject_name, 'pad'),...
%              'signal_start','signal_end','signal_length');
% 
%         fprintf('Completed saving statisics!\n')
%     end
% else    
%    
% end
