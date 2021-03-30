%% Model Script to analyse Anti-Saccade Eye-Link data
% author: thomas.andrillon@monash.edu
clear all
close all

run localdef.m

addpath(genpath(lscpTools_path));
addpath(edf2mat_path);

files=dir([eyet_path filesep '*' filesep '*Acq_clean.mat']);

%% Loop on files
window=[-1000 7000]; %3.4 image alone 2.6 liquid+image
redo=1;
all_CSminus_trials=[];
all_CSplus_trials=[];

all_CSminus_trials_beg=[];
all_CSplus_trials_beg=[];

all_CSminus_trials_end=[];
all_CSplus_trials_end=[];


for n=1:length(files)
    subID=files(n).name(12:15);
    savename=files(n).name;
    
    % load clean data for subID
    load([files(n).folder filesep savename]); %,'EL_headers','EL_data','EL_events');
    
    % find timing of image presentation
    event_types = EL_events.Events.type';
    event_samples = EL_events.Events.time';
    event_CSminus_idx = match_str(event_types,'CSminus');
    event_CSplus_idx = match_str(event_types,'CSplus');
    event_CSminus_samples = event_samples(event_CSminus_idx);
    event_CSplus_samples = event_samples(event_CSplus_idx);
    
    % epoch around image presentation
    CSminus_trials=[];
    for k=1:length(event_CSminus_samples)
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials=[CSminus_trials ; temp_trial];
    end
    
    CSplus_trials=[];
    for k=1:length(event_CSplus_samples)
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials=[CSplus_trials ; temp_trial];
    end
    
    %%%% Just the first 5 trials
    CSminus_trials_beg=[];
    for k=1:5 %length(event_CSminus_samples)
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_beg=[CSminus_trials_beg ; temp_trial];
    end
    
    CSplus_trials_beg=[];
    for k=1:5 %length(event_CSplus_samples)
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_beg=[CSplus_trials_beg ; temp_trial];
    end
    
    %%%% Just the last 5 trials
    CSminus_trials_end=[];
    for k=26:30
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_end=[CSminus_trials_end ; temp_trial];
    end
    
    CSplus_trials_end=[];
    for k=26:30
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_end=[CSplus_trials_end ; temp_trial];
    end
    
    all_CSminus_trials=[all_CSminus_trials ; nanmean(CSminus_trials)];
    all_CSplus_trials=[all_CSplus_trials ; nanmean(CSplus_trials)];
    
    all_CSminus_trials_beg=[all_CSminus_trials_beg ; nanmean(CSminus_trials_beg)];
    all_CSplus_trials_beg=[all_CSplus_trials_beg ; nanmean(CSplus_trials_beg)];
    
    all_CSminus_trials_end=[all_CSminus_trials_end ; nanmean(CSminus_trials_end)];
    all_CSplus_trials_end=[all_CSplus_trials_end ; nanmean(CSplus_trials_end)];
    
end

%%
figure; 
plot(-1000:7000,mean(all_CSminus_trials));
hold on;
plot(-1000:7000,mean(all_CSplus_trials));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
format_fig;
xlim([-1000 7000])

figure; 
subplot(1,2,1);
plot(-1000:7000,mean(all_CSminus_trials_beg));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_beg));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('First 5')
format_fig;
xlim([-1000 7000])

subplot(1,2,2);
plot(-1000:7000,mean(all_CSminus_trials_end));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_end));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('Last 5')
format_fig;
xlim([-1000 7000])
