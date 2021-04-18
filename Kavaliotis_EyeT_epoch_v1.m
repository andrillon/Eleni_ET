%% Model Script to analyse Anti-Saccade Eye-Link data
% author: thomas.andrillon@monash.edu
clear all
close all

run localdef.m

addpath(genpath(lscpTools_path));
addpath(edf2mat_path);

files_A=dir([eyet_path filesep '*' filesep '*Acq_clean.mat']);

%% Loop on files
% Epoching data for the ACQUISITION SESSION

window=[-1000 7000]; %3.4 image alone 2.6 liquid+image
redo=1;

all_CSminus_trials_A=[];
all_CSplus_trials_A=[];

all_CSminus_trials_A1=[];
all_CSplus_trials_A1=[];

all_CSminus_trials_A6=[];
all_CSplus_trials_A6=[];

all_CSminus_trials_A2=[];
all_CSplus_trials_A2=[];

all_CSminus_trials_A3=[];
all_CSplus_trials_A3=[];

all_CSminus_trials_A4=[];
all_CSplus_trials_A4=[];

all_CSminus_trials_A5=[];
all_CSplus_trials_A5=[];


for n=1:length(files_A)
    subID=files_A(n).name(12:15);
    savename=files_A(n).name;
    
    % load clean data for subID
    load([files_A(n).folder filesep savename]); %,'EL_headers','EL_data','EL_events');
    
    % find timing of image presentation
    event_types = EL_events.Events.type';
    event_samples = EL_events.Events.time';
    event_CSminus_idx = match_str(event_types,'CSminus');
    event_CSplus_idx = match_str(event_types,'CSplus');
    event_CSminus_samples = event_samples(event_CSminus_idx);
    event_CSplus_samples = event_samples(event_CSplus_idx);
    
    % Epoch around image presentation for ALL TRIALS
    CSminus_trials_A=[];
    for k=1:length(event_CSminus_samples)
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_A=[CSminus_trials_A ; temp_trial];
    end
    
    CSplus_trials_A=[];
    for k=1:length(event_CSplus_samples)
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_A=[CSplus_trials_A ; temp_trial];
    end
    
    %%%% Trials 1 - 5
    CSminus_trials_A1=[];
    for k=1:5 %length(event_CSminus_samples)
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_A1=[CSminus_trials_A1 ; temp_trial];
    end
    
    CSplus_trials_A1=[];
    for k=1:5 %length(event_CSplus_samples)
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_A1=[CSplus_trials_A1 ; temp_trial];
    end
    
    %%%% Trials 26 - 30
    CSminus_trials_A6=[];
    for k=26:30 %length(event_CSminus_samples)
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_A6=[CSminus_trials_A6 ; temp_trial];
    end
    
    CSplus_trials_A6=[];
    for k=26:30 %length(event_CSplus_samples)
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_A6=[CSplus_trials_A6 ; temp_trial];
    end
    
    %%%% Trials 6 - 10
    
    CSminus_trials_A2=[];
    for k=6:10
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_A2=[CSminus_trials_A2 ; temp_trial];
    end
    
    CSplus_trials_A2=[];
    for k=6:10
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_A2=[CSplus_trials_A2 ; temp_trial];
    end
    
    %%%% Trials 11 - 15
    
    CSminus_trials_A3=[];
    for k=11:15
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_A3=[CSminus_trials_A3 ; temp_trial];
    end
    
    CSplus_trials_A3=[];
    for k=11:15
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_A3=[CSplus_trials_A3 ; temp_trial];
    end
    
     %%%% Trials 16 - 20
    
    CSminus_trials_A4=[];
    for k=16:20
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_A4=[CSminus_trials_A4 ; temp_trial];
    end
    
    CSplus_trials_A4=[];
    for k=16:20
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_A4=[CSplus_trials_A4 ; temp_trial];
    end
      %%%% Trials 21 - 25
    
    CSminus_trials_A5=[];
    for k=21:25
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_A5=[CSminus_trials_A5 ; temp_trial];
    end
    
    CSplus_trials_A5=[];
    for k=21:25
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_A5=[CSplus_trials_A5 ; temp_trial];
    end
         
    all_CSminus_trials_A=[all_CSminus_trials_A ; nanmean(CSminus_trials_A)];
    all_CSplus_trials_A=[all_CSplus_trials_A ; nanmean(CSplus_trials_A)];
    
    all_CSminus_trials_A1=[all_CSminus_trials_A1 ; nanmean(CSminus_trials_A1)];
    all_CSplus_trials_A1=[all_CSplus_trials_A1 ; nanmean(CSplus_trials_A1)];
    
    all_CSminus_trials_A6=[all_CSminus_trials_A6 ; nanmean(CSminus_trials_A6)];
    all_CSplus_trials_A6=[all_CSplus_trials_A6 ; nanmean(CSplus_trials_A6)];
    
    all_CSminus_trials_A2=[all_CSminus_trials_A2 ; nanmean(CSminus_trials_A2)];
    all_CSplus_trials_A2=[all_CSplus_trials_A2 ; nanmean(CSplus_trials_A2)];
    
    all_CSminus_trials_A3=[all_CSminus_trials_A3 ; nanmean(CSminus_trials_A3)];
    all_CSplus_trials_A3=[all_CSplus_trials_A3 ; nanmean(CSplus_trials_A3)];
    
    all_CSminus_trials_A4=[all_CSminus_trials_A4 ; nanmean(CSminus_trials_A4)];
    all_CSplus_trials_A4=[all_CSplus_trials_A4 ; nanmean(CSplus_trials_A4)];
    
    all_CSminus_trials_A5=[all_CSminus_trials_A5 ; nanmean(CSminus_trials_A5)];
    all_CSplus_trials_A5=[all_CSplus_trials_A5 ; nanmean(CSplus_trials_A5)];
end

%%
%%% Average of ALL trials
figure(8); 
plot(-1000:7000,mean(all_CSminus_trials_A));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_A));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
format_fig;
xlim([-1000 7000])
%%
%%% Comparing first and last 5 trials 
figure(7); 
subplot(1,2,1);
plot(-1000:7000,mean(all_CSminus_trials_A1));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_A1));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('A1 - Trials 1-5')
format_fig;
xlim([-1000 7000])

subplot(1,2,2);
plot(-1000:7000,mean(all_CSminus_trials_A6));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_A6));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('A6 - Trials 26-30')
format_fig;
xlim([-1000 7000])
%%
%%% Average of each section of trials

figure (1); 
plot(-1000:7000,mean(all_CSminus_trials_A1));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_A1));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('A1 - Trials 1-5')
format_fig;
xlim([-1000 7000])
ylim([0 350])

figure (2);
plot(-1000:7000,mean(all_CSminus_trials_A2));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_A2));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('A2 - Trials 6-10')
format_fig;
xlim([-1000 7000])
ylim([0 350])

figure (3); 
plot(-1000:7000,mean(all_CSminus_trials_A3));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_A3));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('A3 - Trials 11-15')
format_fig;
xlim([-1000 7000])
ylim([0 350])

figure (4);
plot(-1000:7000,mean(all_CSminus_trials_A4));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_A4));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('A4 - Trials 16-20')
format_fig;
xlim([-1000 7000])
ylim([0 350])

figure (5);
plot(-1000:7000,mean(all_CSminus_trials_A5));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_A5));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('A5 - Trials 21-25')
format_fig;
xlim([-1000 7000])
ylim([0 350])

figure(6);
plot(-1000:7000,mean(all_CSminus_trials_A6));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_A6));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('A6 - Trials 26-30')
format_fig;
xlim([-1000 7000])
ylim([0 350])

%%
files_E=dir([eyet_path filesep '*' filesep '*Ext_clean.mat']);

%% Loop on files
% Epoching data for the EXTINCTION SESSION

window=[-1000 7000]; %3.4 image alone 2.6 liquid+image
redo=1;

all_CSminus_trials_E=[];
all_CSplus_trials_E=[];

all_CSminus_trials_E1=[];
all_CSplus_trials_E1=[];

all_CSminus_trials_E6=[];
all_CSplus_trials_E6=[];

all_CSminus_trials_E2=[];
all_CSplus_trials_E2=[];

all_CSminus_trials_E3=[];
all_CSplus_trials_E3=[];

all_CSminus_trials_E4=[];
all_CSplus_trials_E4=[];

all_CSminus_trials_E5=[];
all_CSplus_trials_E5=[];


for n=1:length(files_E)
    subID=files_E(n).name(12:15);
    savename=files_E(n).name;
    
    % load clean data for subID
    load([files_E(n).folder filesep savename]); %,'EL_headers','EL_data','EL_events');
    
    % find timing of image presentation
    event_types = EL_events.Events.type';
    event_samples = EL_events.Events.time';
    event_CSminus_idx = match_str(event_types,'CSminus');
    event_CSplus_idx = match_str(event_types,'CSplus');
    event_CSminus_samples = event_samples(event_CSminus_idx);
    event_CSplus_samples = event_samples(event_CSplus_idx);
    
    % Epoch around image presentation for ALL TRIALS
    CSminus_trials_E=[];
    for k=1:length(event_CSminus_samples)
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_E=[CSminus_trials_E; temp_trial];
    end
    
    CSplus_trials_E=[];
    for k=1:length(event_CSplus_samples)
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_E=[CSplus_trials_E ; temp_trial];
    end
    
    %%%% Trials 1 - 5
    CSminus_trials_E1=[];
    for k=1:5 %length(event_CSminus_samples)
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_E1=[CSminus_trials_E1 ; temp_trial];
    end
    
    CSplus_trials_E1=[];
    for k=1:5 %length(event_CSplus_samples)
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_E1=[CSplus_trials_E1 ; temp_trial];
    end
    
    %%%% Trials 26 - 30
    CSminus_trials_E6=[];
    for k=26:30 %length(event_CSminus_samples)
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_E6=[CSminus_trials_E6 ; temp_trial];
    end
    
    CSplus_trials_E6=[];
    for k=26:30 %length(event_CSplus_samples)
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_E6=[CSplus_trials_E6 ; temp_trial];
    end
    
    %%%% Trials 6 - 10
    
    CSminus_trials_E2=[];
    for k=6:10
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_E2=[CSminus_trials_E2 ; temp_trial];
    end
    
    CSplus_trials_E2=[];
    for k=6:10
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_E2=[CSplus_trials_E2 ; temp_trial];
    end
    
    %%%% Trials 11 - 15
    
    CSminus_trials_E3=[];
    for k=11:15
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_E3=[CSminus_trials_E3 ; temp_trial];
    end
    
    CSplus_trials_E3=[];
    for k=11:15
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_E3=[CSplus_trials_E3 ; temp_trial];
    end
    
     %%%% Trials 16 - 20
    
    CSminus_trials_E4=[];
    for k=16:20
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_E4=[CSminus_trials_E4 ; temp_trial];
    end
    
    CSplus_trials_E4=[];
    for k=16:20
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_E4=[CSplus_trials_E4 ; temp_trial];
    end
      %%%% Trials 21 - 25
    
    CSminus_trials_E5=[];
    for k=21:25
        temp_trialonset = find(EL_data.time==event_CSminus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-mean(temp_trial(1:1000));
        CSminus_trials_E5=[CSminus_trials_E5 ; temp_trial];
    end
    
    CSplus_trials_E5=[];
    for k=21:25
        temp_trialonset = find(EL_data.time==event_CSplus_samples(k));
        temp_trial=EL_data.clean_pupilSize((temp_trialonset+window(1)):(temp_trialonset+window(2)),1)';
        % substract the baseline
        temp_trial=temp_trial-nanmean(temp_trial(1:1000));
        CSplus_trials_E5=[CSplus_trials_E5 ; temp_trial];
    end
         
    all_CSminus_trials_E=[all_CSminus_trials_E ; nanmean(CSminus_trials_E)];
    all_CSplus_trials_E=[all_CSplus_trials_E ; nanmean(CSplus_trials_E)];
    
    all_CSminus_trials_E1=[all_CSminus_trials_E1 ; nanmean(CSminus_trials_E1)];
    all_CSplus_trials_E1=[all_CSplus_trials_E1 ; nanmean(CSplus_trials_E1)];
    
    all_CSminus_trials_E6=[all_CSminus_trials_E6 ; nanmean(CSminus_trials_E6)];
    all_CSplus_trials_E6=[all_CSplus_trials_E6 ; nanmean(CSplus_trials_E6)];
    
    all_CSminus_trials_E2=[all_CSminus_trials_E2 ; nanmean(CSminus_trials_E2)];
    all_CSplus_trials_E2=[all_CSplus_trials_E2 ; nanmean(CSplus_trials_E2)];
    
    all_CSminus_trials_E3=[all_CSminus_trials_E3 ; nanmean(CSminus_trials_E3)];
    all_CSplus_trials_E3=[all_CSplus_trials_E3 ; nanmean(CSplus_trials_E3)];
    
    all_CSminus_trials_E4=[all_CSminus_trials_E4 ; nanmean(CSminus_trials_E4)];
    all_CSplus_trials_E4=[all_CSplus_trials_E4 ; nanmean(CSplus_trials_E4)];
    
    all_CSminus_trials_E5=[all_CSminus_trials_E5 ; nanmean(CSminus_trials_E5)];
    all_CSplus_trials_E5=[all_CSplus_trials_E5 ; nanmean(CSplus_trials_E5)];
end

%%
%%% Average of ALL trials
figure(18); 
plot(-1000:7000,mean(all_CSminus_trials_E));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_E));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
format_fig;
xlim([-1000 7000])
%%
%%% Comparing first and last 5 trials 
figure(17); 
subplot(1,2,1);
plot(-1000:7000,mean(all_CSminus_trials_E1));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_E1));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('E1 - Trials 1-5')
format_fig;
xlim([-1000 7000])

subplot(1,2,2);
plot(-1000:7000,mean(all_CSminus_trials_E6));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_E6));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('E6 - Trials 26-30')
format_fig;
xlim([-1000 7000])
%%
%%% Average of each section of trials

figure (11); 
plot(-1000:7000,mean(all_CSminus_trials_E1));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_E1));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('E1 - Trials 1-5')
format_fig;
xlim([-1000 7000])
ylim([0 350])

figure (12);
plot(-1000:7000,mean(all_CSminus_trials_E2));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_E2));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('E2 - Trials 6-10')
format_fig;
xlim([-1000 7000])
ylim([0 350])

figure (13); 
plot(-1000:7000,mean(all_CSminus_trials_E3));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_E3));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('E3 - Trials 11-15')
format_fig;
xlim([-1000 7000])
ylim([0 350])

figure (14);
plot(-1000:7000,mean(all_CSminus_trials_E4));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_E4));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('E4 - Trials 16-20')
format_fig;
xlim([-1000 7000])
ylim([0 350])

figure (15);
plot(-1000:7000,mean(all_CSminus_trials_E5));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_E5));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('E5 - Trials 21-25')
format_fig;
xlim([-1000 7000])
ylim([0 350])

figure(16);
plot(-1000:7000,mean(all_CSminus_trials_E6));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_E6));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('E6 - Trials 26-30')
format_fig;
xlim([-1000 7000])
ylim([0 350])
%%


