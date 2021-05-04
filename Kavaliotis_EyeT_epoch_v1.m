%% Model Script to analyse Anti-Saccade Eye-Link data
% author: thomas.andrillon@monash.edu

%SECTION 1
clear all
close all

run localdef.m

addpath(genpath(lscpTools_path));
addpath(edf2mat_path);

files_A=dir([eyet_path filesep '*' filesep '*Acq_clean.mat']);

%% Loop on files
% SECTION 2: Epoching data for the ACQUISITION SESSION

window=[-1000 7000]; %2.4 image alone 3.6 liquid+image
windowtimes=(window(1):window(2))/1000;
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

all_Names=[];

all_CSminus_avValues=[];
all_CSplus_avValues=[];

all_CSminus_avValues_A=[];
all_CSplus_avValues_A=[];

for n=1:length(files_A)
    subID=files_A(n).name(12:15);
    savename=files_A(n).name;
    all_Names=[all_Names {subID}];
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
    
    all_CSminus_avValues=[all_CSminus_avValues ; [str2num(subID(2:end))*ones(size(CSminus_trials_A,1),1) (1:size(CSminus_trials_A,1))' [max(CSminus_trials_A(:,windowtimes>0 & windowtimes<2.4),[],2)  median(CSminus_trials_A(:,windowtimes>0 & windowtimes<2.4),2) mean(CSminus_trials_A(:,windowtimes>0 & windowtimes<2.4),2) ...
        max(CSminus_trials_A(:,windowtimes>2.4 & windowtimes<2.4+3.6),[],2)  median(CSminus_trials_A(:,windowtimes>2.4 & windowtimes<2.4+3.6),2) mean(CSminus_trials_A(:,windowtimes>2.4 & windowtimes<2.4+3.6),2)]]];

    all_CSplus_avValues=[all_CSplus_avValues ; [str2num(subID(2:end))*ones(size(CSplus_trials_A,1),1) (1:size(CSplus_trials_A,1))'  [max(CSplus_trials_A(:,windowtimes>0 & windowtimes<2.4),[],2)  median(CSplus_trials_A(:,windowtimes>0 & windowtimes<2.4),2) mean(CSplus_trials_A(:,windowtimes>0 & windowtimes<2.4),2) ...
        max(CSplus_trials_A(:,windowtimes>2.4 & windowtimes<2.4+3.6),[],2)  median(CSplus_trials_A(:,windowtimes>2.4 & windowtimes<2.4+3.6),2) mean(CSplus_trials_A(:,windowtimes>3.4 & windowtimes<2.4+3.6),2)]]];    
    
    % This avValues is not split between image alone and image+liquid, it
    % is entire CS presentation.
    
    all_CSminus_avValues_A=[all_CSminus_avValues_A ; [str2num(subID(2:end))*ones(size(CSminus_trials_A,1),1) (1:size(CSminus_trials_A,1))' [max(CSminus_trials_A(:,windowtimes>0 & windowtimes<6.0),[],2)  median(CSminus_trials_A(:,windowtimes>0 & windowtimes<6.0),2) mean(CSminus_trials_A(:,windowtimes>0 & windowtimes<6.0),2)]]];
    
    all_CSplus_avValues_A=[all_CSplus_avValues_A ; [str2num(subID(2:end))*ones(size(CSplus_trials_A,1),1) (1:size(CSplus_trials_A,1))'  [max(CSplus_trials_A(:,windowtimes>0 & windowtimes<6.0),[],2)  median(CSplus_trials_A(:,windowtimes>0 & windowtimes<6.0),2) mean(CSplus_trials_A(:,windowtimes>0 & windowtimes<6.0),2)]]];    

end

%% SECTION 3: Create a table for the Mean, Median and Mode for CSplus and CSminus

for n=1:length(files_A)

    all_CSminus_avValues_table=array2table(all_CSminus_avValues,'VariableNames',{'SubID','TrialNumber','Max_WindowImage','Median_WindowImage','Mean_WindowImage',...
    'Max_WindowImageLiquid','Median_WindowImageLiquid','Mean_WindowImageLiquid'});

    all_CSplus_avValues_table=array2table(all_CSplus_avValues,'VariableNames',{'SubID','TrialNumber','Max_WindowImage','Median_WindowImage','Mean_WindowImage',...
    'Max_WindowImageLiquid','Median_WindowImageLiquid','Mean_WindowImageLiquid'});

    all_CSminus_avValues_A_table=array2table(all_CSminus_avValues_A,'VariableNames',{'SubID','TrialNumber','Max_Window','Median_Window','Mean_Window'});

    all_CSplus_avValues_A_table=array2table(all_CSplus_avValues_A,'VariableNames',{'SubID','TrialNumber','Max_Window','Median_Window','Mean_Window'});
end 

%% SECTION 4: Plotting average pupil size 

%%% Plot average of ALL trials in the same graph
figure(401); 
plot(-1000:7000,mean(all_CSminus_trials_A));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_A));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
format_fig;
xlim([-1000 7000])
ylim([0 300])

%%%% Plot average of all trials + the SEM for CSplus and minus in the same
%%%% and seperate graphs
figure(402); 
simpleTplot(-1000:7000,all_CSminus_trials_A,0,'b',0,'-',0.5,1,0,1,3);
legend({'SEM','CS-'})
xlabel('Time (s)')
ylabel('Pupil size')
title('CS-')
format_fig;
xlim([-1000 7000])

figure(403); 
simpleTplot(-1000:7000,all_CSplus_trials_A,0,'r',0,'-',0.5,1,0,1,3);
legend({'SEM','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('CS+')
format_fig;
xlim([-1000 7000])

figure (404);
simpleTplot(-1000:7000,all_CSminus_trials_A,0,'b',0,'-',0.5,1,0,1,3);
hold on ;
simpleTplot(-1000:7000,all_CSplus_trials_A,0,'r',0,'-',0.5,1,0,1,3);
legend({'CS- SEM','CS-','CS+ SEM', 'CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('Mean and SEM for CS+ and CS-')
format_fig;
xlim([-1000 7000])

%%% Plot average of ALL trials for EACH participant in the same graph for
%%% CS- and CS+ seperately

figure(405); format_fig;
plot(-1000:7000,all_CSminus_trials_A');
title('CS- Trials')
legend(all_Names);
xlim([-1000 7000])
ylim([-100 400])

figure(406); format_fig;
plot(-1000:7000,all_CSplus_trials_A');
title('CS+ Trials')
legend(all_Names);
xlim([-1000 7000])
ylim([-100 400])

%% SECTION 5: Average across subjects for trial numbers

uniqueSubIDs=unique(all_CSminus_avValues(:,1));
uniqueSubIDs_A=unique(all_CSminus_avValues_A(:,1));

%Splitting mean pupil size at each trial between image alone and with liquid and by
%participant for CS-:

for k=1:length(uniqueSubIDs)
figure(501);
subplot(1,length(uniqueSubIDs),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),5),all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),2))) %%Mean value for image alone in CS- trials
plot(1:30,grpstats(all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),8),all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),2))) %%Mean value for image+liquid in CS- trials 
title(k)
legend({'CS-alone','CS-with liquid'})
ylim([0 600])
end

%Splitting mean pupil size at each trial between image alone and with liquid and by
%participant for CS+:

for k=1:length(uniqueSubIDs)
figure(502);
subplot(1,length(uniqueSubIDs),k);
hold on;
plot(1:30,grpstats(all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),5),all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),2))) %%Mean value for image alone in CS+ trials
plot(1:30,grpstats(all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),8),all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),2))) %%Mean value for image+liquid in CS+ trials 
title(k)
legend({'CS+alone','CS+with liquid'})
ylim([0 600])
end

%Splitting mean pupil size at each trial by each CS type and by participant
%for entire image presentation (e.g., alone and with liquid):
for k=1:length(uniqueSubIDs_A)
figure(503);
subplot(1,length(uniqueSubIDs_A),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),5),all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Mean value (e.g., column 5) for entire CS- presentation
hold on;
plot(1:30,grpstats(all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),5),all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Mean value (e.g., column 5) for entire CS+ presentation
title(k)
ylim([0 600])
end

%With more detail to graph, same as above:
for k=1:length(uniqueSubIDs_A)
figure (504);
subplot(1,length(uniqueSubIDs_A),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),5),all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Mean value (e.g., cloumn 5) for entire CS- presentation
hold on;
plot(1:30,grpstats(all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),5),all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Mean value (e.g., cloumn 5) for entire CS+ presentation
title(k)
legend({'CS-','CS+'})
xlabel('Trials')
ylabel('Pupil size')
format_fig;
ylim([0 600])
end

%Splitting max pupil size at each trial by each CS type and by participant
%for entire image presentation (e.g., alone and with liquid):

for k=1:length(uniqueSubIDs_A)
figure (505);
subplot(1,length(uniqueSubIDs_A),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),3),all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Max value (e.g., cloumn 5) for entire CS- presentation
hold on;
plot(1:30,grpstats(all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),3),all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Max value (e.g., cloumn 5) for entire CS+ presentation
title(k)
legend({'CS-','CS+'})
xlabel('Trials')
ylabel('Pupil size')
format_fig;
end

%Splitting median pupil size at each trial by each CS type and by participant
%for entire image presentation (e.g., alone and with liquid):

for k=1:length(uniqueSubIDs_A)
figure (506);
subplot(1,length(uniqueSubIDs_A),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),4),all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Median value (e.g., cloumn 5) for entire CS- presentation
hold on;
plot(1:30,grpstats(all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),4),all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Median value (e.g., cloumn 5) for entire CS+ presentation
title(k)
legend({'CS-','CS+'})
xlabel('Trials')
ylabel('Pupil size')
format_fig;
end

%Split out CS+ and CS- for max scores
for k=1:length(uniqueSubIDs_A)
figure (507);
subplot(1,length(uniqueSubIDs_A),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),3),all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Mean value (e.g., cloumn 5) for entire CS- presentation
title(k)
sgtitle('CS-','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;
end

for k=1:length(uniqueSubIDs_A)
figure (508);
subplot(1,length(uniqueSubIDs_A),k);
hold on;
plot(1:30,grpstats(all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),3),all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Mean value (e.g., cloumn 5) for entire CS+ presentation
title(k)
sgtitle('CS+','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;
end

%TRIAL BY TRIAL FOR OVERALL SCORES FOR MEDIAN, MAX AND MEAN for entire image presentation (e.g., alone and with liquid):
% MEAN
figure (509);
plot(1:30,grpstats(all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),5),all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Mean value (e.g., cloumn 5) for entire CS- presentation
title('CS-','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

figure (510);
plot(1:30,grpstats(all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),5),all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Mean value (e.g., cloumn 5) for entire CS+ presentation
title('CS+','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

figure (5091);
plot(1:30,grpstats(all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),5),all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Mean value (e.g., cloumn 5) for entire CS- presentation
hold on
plot(1:30,grpstats(all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),5),all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Mean value (e.g., cloumn 5) for entire CS+ presentation
xlabel('Trials')
ylabel('Pupil size')
legend('CS-', 'CS+')
format_fig;

% MEDIAN
figure (511);
plot(1:30,grpstats(all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),4),all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Median value (e.g., cloumn 4) for entire CS- presentation
title('CS-','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

figure (512);
plot(1:30,grpstats(all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),4),all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Median value (e.g., cloumn 4) for entire CS+ presentation
title('CS+','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

%MAX
figure (513);
plot(1:30,grpstats(all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),3),all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Max value (e.g., cloumn 3) for entire CS- presentation
title('CS-','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

figure (514);
plot(1:30,grpstats(all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),3),all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Max value (e.g., cloumn 3) for entire CS+ presentation
title('CS+','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

figure (5141);
plot(1:30,grpstats(all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),3),all_CSminus_avValues_A(all_CSminus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Max value (e.g., cloumn 3) for entire CS- presentation
hold on
plot(1:30,grpstats(all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),3),all_CSplus_avValues_A(all_CSplus_avValues_A(:,1)==uniqueSubIDs_A(k),2)))% Max value (e.g., cloumn 3) for entire CS+ presentation
legend ('CS-', 'CS+')
xlabel('Trials')
ylabel('Pupil size')
format_fig;


%TRIAL BY TRIAL FOR OVERALL SCORES FOR MEDIAN, MAX AND MEAN for split of
%image alone and with liquiD:
%MEAN 
figure(515);
plot(1:30,grpstats(all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),5),all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),2)))
hold on
plot(1:30,grpstats(all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),8),all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),2))) 
legend({'CS-alone','CS-image+liquid'})

figure (516);
plot(1:30,grpstats(all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),5),all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),2))) 
hold on
plot(1:30,grpstats(all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),8),all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),2)))  
legend({'CS+alone','CS+image+liquid'})

%MEDIAN for split of image alone and with liquid
figure(517);
plot(1:30,grpstats(all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),4),all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),2))) 
hold on
plot(1:30,grpstats(all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),7),all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),2))) 

figure (518);
plot(1:30,grpstats(all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),4),all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),2)))
hold on
plot(1:30,grpstats(all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),7),all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),2)))  

%MAX for split of image alone and with liquid
figure(519);
plot(1:30,grpstats(all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),3),all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),2))) 
hold on
plot(1:30,grpstats(all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),6),all_CSminus_avValues(all_CSminus_avValues(:,1)==uniqueSubIDs(k),2))) 
legend({'CS-alone','CS-image+liquid'})

figure (520);
plot(1:30,grpstats(all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),3),all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),2))) 
hold on
plot(1:30,grpstats(all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),6),all_CSplus_avValues(all_CSplus_avValues(:,1)==uniqueSubIDs(k),2))) 
legend({'CS+alone','CS+image+liquid'})

%%
%%% SECTION 6: Comparing first and last 5 trials 
figure(601); 
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
ylim([0 400])

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
ylim([0 400])

%%
%%% SECTION 7: Average of each section of trials

%A1: 
figure (701); 
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

%For each pt average of A1 for each CS type on seperate graphs:

figure(7011); format_fig;
plot(-1000:7000,all_CSminus_trials_A1');
title('A1 - Trials 1-5 for CS-')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

figure(7012); format_fig;
plot(-1000:7000,all_CSplus_trials_A1');
title('A1 - Trials 1-5 for CS+')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

%A2

figure (702);
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

figure(7021); format_fig;
plot(-1000:7000,all_CSminus_trials_A2');
title('A2 - Trials 6-10 for CS-')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

figure(7022); format_fig;
plot(-1000:7000,all_CSplus_trials_A2');
title('A2 - Trials 6-10 for CS+')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

%A3:

figure (703); 
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

figure(7031); format_fig;
plot(-1000:7000,all_CSminus_trials_A3');
title('A3 - Trials 11-15 for CS-')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

figure(7032); format_fig;
plot(-1000:7000,all_CSplus_trials_A3');
title('A3 - Trials 11-15 for CS+')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

%A4:

figure (704);
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

figure(7041); format_fig;
plot(-1000:7000,all_CSminus_trials_A4');
title('A4 - Trials 16-20 for CS-')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

figure(7042); format_fig;
plot(-1000:7000,all_CSplus_trials_A4');
title('A4 - Trials 16-20 for CS+')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

%A5:

figure (705);
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

figure(7051); format_fig;
plot(-1000:7000,all_CSminus_trials_A5');
title('A5 - Trials 20-25 for CS-')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

figure(7052); format_fig;
plot(-1000:7000,all_CSplus_trials_A5');
title('A5 - Trials 20-25 for CS+')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

%A6: 

figure(706);
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

figure(7061); format_fig;
plot(-1000:7000,all_CSminus_trials_A6');
title('A4 - Trials 26-30 for CS-')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

figure(7062); format_fig;
plot(-1000:7000,all_CSplus_trials_A6');
title('A6 - Trials 26-30 for CS+')
legend(all_Names);
xlim([-1000 7000])
ylim([0 600])

%% SECTION 8: EXTINCTION DATA:

files_E=dir([eyet_path filesep '*' filesep '*Ext_clean.mat']);

%% Loop on files
% SECTION 9: Epoching data for the EXTINCTION SESSION

window=[-1000 7000]; %2.4 image alone 3.6 liquid+image
windowtimes=(window(1):window(2))/1000;
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

all_Names=[];

all_CSminus_avValuesE=[];
all_CSplus_avValuesE=[];

all_CSminus_avValues_Ext=[];
all_CSplus_avValues_Ext=[];



for n=1:length(files_E)
    subID=files_E(n).name(12:15);
    savename=files_E(n).name;
    all_Names=[all_Names {subID}];
    
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
    
    all_CSminus_avValuesE=[all_CSminus_avValuesE ; [str2num(subID(2:end))*ones(size(CSminus_trials_E,1),1) (1:size(CSminus_trials_E,1))' [max(CSminus_trials_E(:,windowtimes>0 & windowtimes<2.4),[],2)  median(CSminus_trials_E(:,windowtimes>0 & windowtimes<2.4),2) mean(CSminus_trials_E(:,windowtimes>0 & windowtimes<2.4),2) ...
        max(CSminus_trials_E(:,windowtimes>2.4 & windowtimes<2.4+3.6),[],2)  median(CSminus_trials_E(:,windowtimes>2.4 & windowtimes<2.4+3.6),2) mean(CSminus_trials_E(:,windowtimes>2.4 & windowtimes<2.4+3.6),2)]]];

    all_CSplus_avValuesE=[all_CSplus_avValuesE ; [str2num(subID(2:end))*ones(size(CSplus_trials_E,1),1) (1:size(CSplus_trials_E,1))'  [max(CSplus_trials_E(:,windowtimes>0 & windowtimes<2.4),[],2)  median(CSplus_trials_E(:,windowtimes>0 & windowtimes<2.4),2) mean(CSplus_trials_E(:,windowtimes>0 & windowtimes<2.4),2) ...
        max(CSplus_trials_E(:,windowtimes>2.4 & windowtimes<2.4+3.6),[],2)  median(CSplus_trials_E(:,windowtimes>2.4 & windowtimes<2.4+3.6),2) mean(CSplus_trials_E(:,windowtimes>3.4 & windowtimes<2.4+3.6),2)]]];    
    
    % This avValues is not split between image alone and image+liquid, it
    % is entire CS presentation.
    
    all_CSminus_avValues_Ext=[all_CSminus_avValues_Ext ; [str2num(subID(2:end))*ones(size(CSminus_trials_E,1),1) (1:size(CSminus_trials_E,1))' [max(CSminus_trials_E(:,windowtimes>0 & windowtimes<6.0),[],2)  median(CSminus_trials_E(:,windowtimes>0 & windowtimes<6.0),2) mean(CSminus_trials_E(:,windowtimes>0 & windowtimes<6.0),2)]]];
    
    all_CSplus_avValues_Ext=[all_CSplus_avValues_Ext ; [str2num(subID(2:end))*ones(size(CSplus_trials_E,1),1) (1:size(CSplus_trials_E,1))'  [max(CSplus_trials_E(:,windowtimes>0 & windowtimes<6.0),[],2)  median(CSplus_trials_E(:,windowtimes>0 & windowtimes<6.0),2) mean(CSplus_trials_E(:,windowtimes>0 & windowtimes<6.0),2)]]];    

end

%% SECTION 10: Create a table for the Mean, Median and Mode for CSplus and CSminus for EXTINCTION TRIALS
for n=1:length(files_E)

    all_CSminus_avValuesE_table=array2table(all_CSminus_avValuesE,'VariableNames',{'SubID','TrialNumber','Max_WindowImage','Median_WindowImage','Mean_WindowImage',...
    'Max_WindowImageLiquid','Median_WindowImageLiquid','Mean_WindowImageLiquid'});

    all_CSplus_avValuesE_table=array2table(all_CSplus_avValuesE,'VariableNames',{'SubID','TrialNumber','Max_WindowImage','Median_WindowImage','Mean_WindowImage',...
    'Max_WindowImageLiquid','Median_WindowImageLiquid','Mean_WindowImageLiquid'});

    all_CSminus_avValues_Ext_table=array2table(all_CSminus_avValues_Ext,'VariableNames',{'SubID','TrialNumber','Max_Window','Median_Window','Mean_Window'});

    all_CSplus_avValues_Ext_table=array2table(all_CSplus_avValues_Ext,'VariableNames',{'SubID','TrialNumber','Max_Window','Median_Window','Mean_Window'});
end 


%%
%%% SECTION 11: Average of ALL trials
figure(1101); 
plot(-1000:7000,mean(all_CSminus_trials_E));
hold on;
plot(-1000:7000,mean(all_CSplus_trials_E));
legend({'CS-','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
format_fig;
xlim([-1000 7000])


%%%% Plot average of all trials + the SEM for CSplus and minus in the same
%%%% and seperate graphs
figure(1102); 
simpleTplot(-1000:7000,all_CSminus_trials_E,0,'b',0,'-',0.5,1,0,1,3);
legend({'SEM','CS-'})
xlabel('Time (s)')
ylabel('Pupil size')
title('CS-')
format_fig;
xlim([-1000 7000])

figure(1103); 
simpleTplot(-1000:7000,all_CSplus_trials_E,0,'r',0,'-',0.5,1,0,1,3);
legend({'SEM','CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('CS+')
format_fig;
xlim([-1000 7000])

figure (1104);
simpleTplot(-1000:7000,all_CSminus_trials_E,0,'b',0,'-',0.5,1,0,1,3);
hold on ;
simpleTplot(-1000:7000,all_CSplus_trials_E,0,'r',0,'-',0.5,1,0,1,3);
legend({'CS- SEM','CS-','CS+ SEM', 'CS+'})
xlabel('Time (s)')
ylabel('Pupil size')
title('Mean and SEM for CS+ and CS-')
format_fig;
xlim([-1000 7000])

%%% Plot average of ALL trials for each participant in the same graph for
%%% CS- and CS+ seperately

figure(1105); format_fig;
plot(-1000:7000,all_CSminus_trials_E');
title('CS- Trials')
legend(all_Names, 'Location','northwest');
ylim([-200 350])

figure(1106); format_fig;
plot(-1000:7000,all_CSplus_trials_E');
title('CS+ Trials')
legend(all_Names, 'Location','northwest');
ylim([-200 350])



%% SECTION 12: Average across subjects for trial numbers

uniqueSubIDsE=unique(all_CSminus_avValuesE(:,1));
uniqueSubIDs_Ext=unique(all_CSminus_avValues_Ext(:,1));

%Splitting mean pupil size at each trial between image alone and with liquid and by
%participant for CS-:

for k=1:length(uniqueSubIDsE)
figure(1201);
subplot(1,length(uniqueSubIDsE),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),5),all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),2))) %%Mean value for image alone in CS- trials
plot(1:30,grpstats(all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),8),all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),2))) %%Mean value for image+liquid in CS- trials 
title(k)
legend({'CS-alone','CS-with liquid'})
ylim([-800 800])
end

%Splitting mean pupil size at each trial between image alone and with liquid and by
%participant for CS+:

for k=1:length(uniqueSubIDsE)
figure(1202);
subplot(1,length(uniqueSubIDsE),k);
hold on;
plot(1:30,grpstats(all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),5),all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),2))) %%Mean value for image alone in CS+ trials
plot(1:30,grpstats(all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),8),all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),2))) %%Mean value for image+liquid in CS+ trials 
title(k)
legend({'CS+alone','CS+with liquid'})
ylim([-500 800])
end

%Splitting mean pupil size at each trial by each CS type and by participant
%for entire image presentation (e.g., alone and with liquid):
for k=1:length(uniqueSubIDs_Ext)
figure(1203);
subplot(1,length(uniqueSubIDs_Ext),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),5),all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Mean value (e.g., column 5) for entire CS- presentation
hold on;
plot(1:30,grpstats(all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),5),all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Mean value (e.g., column 5) for entire CS+ presentation
title(k)
end

%With more detail to graph, same as above:
for k=1:length(uniqueSubIDs_Ext)
figure (1204);
subplot(1,length(uniqueSubIDs_Ext),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),5),all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Mean value (e.g., column 5) for entire CS- presentation
hold on;
plot(1:30,grpstats(all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),5),all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Mean value (e.g., column 5) for entire CS+ presentation
title(k)
legend({'CS-','CS+'})
xlabel('Trials')
ylabel('Pupil size')
ylim([-400 800])
format_fig;
end

%Splitting max pupil size at each trial by each CS type and by participant
%for entire image presentation (e.g., alone and with liquid):

for k=1:length(uniqueSubIDs_Ext)
figure (1205);
subplot(1,length(uniqueSubIDs_Ext),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),3),all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Max value (e.g., cloumn 5) for entire CS- presentation
hold on;
plot(1:30,grpstats(all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),3),all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Max value (e.g., cloumn 5) for entire CS+ presentation
title(k)
legend({'CS-','CS+'})
xlabel('Trials')
ylabel('Pupil size')
format_fig;
end

%Splitting median pupil size at each trial by each CS type and by participant
%for entire image presentation (e.g., alone and with liquid):

for k=1:length(uniqueSubIDs_Ext)
figure (1206);
subplot(1,length(uniqueSubIDs_Ext),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),4),all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Median value (e.g., cloumn 5) for entire CS- presentation
hold on;
plot(1:30,grpstats(all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),4),all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Median value (e.g., cloumn 5) for entire CS+ presentation
title(k)
legend({'CS-','CS+'})
xlabel('Trials')
ylabel('Pupil size')
format_fig;
end

%Split out CS+ and CS- for max scores
for k=1:length(uniqueSubIDs_Ext)
figure (1207);
subplot(1,length(uniqueSubIDs_Ext),k);
hold on;
plot(1:30,grpstats(all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),3),all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Mean value (e.g., cloumn 5) for entire CS- presentation
title(k)
sgtitle('CS-','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;
end

for k=1:length(uniqueSubIDs_Ext)
figure (1208);
subplot(1,length(uniqueSubIDs_Ext),k);
hold on;
plot(1:30,grpstats(all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),3),all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Mean value (e.g., cloumn 5) for entire CS+ presentation
title(k)
sgtitle('CS+','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;
end

%TRIAL BY TRIAL FOR OVERALL SCORES FOR MEDIAN, MAX AND MEAN for entire image presentation (e.g., alone and with liquid):
% MEAN
figure (1209);
plot(1:30,grpstats(all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),5),all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Mean value (e.g., cloumn 5) for entire CS- presentation
title('CS-','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

figure (1210);
plot(1:30,grpstats(all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),5),all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Mean value (e.g., cloumn 5) for entire CS+ presentation
title('CS+','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

figure (12091);
plot(1:30,grpstats(all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),5),all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Mean value (e.g., cloumn 5) for entire CS- presentation
hold on
plot(1:30,grpstats(all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),5),all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Mean value (e.g., cloumn 5) for entire CS+ presentation
legend('CS-','CS+')
xlabel('Trials')
ylabel('Pupil size')
format_fig;

% MEDIAN
figure (1211);
plot(1:30,grpstats(all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),4),all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Median value (e.g., cloumn 4) for entire CS- presentation
title('CS-','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

figure (1212);
plot(1:30,grpstats(all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),4),all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Median value (e.g., cloumn 4) for entire CS+ presentation
title('CS+','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

%MAX
figure (1213);
plot(1:30,grpstats(all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),3),all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Max value (e.g., cloumn 3) for entire CS- presentation
title('CS-','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

figure (1214);
plot(1:30,grpstats(all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),3),all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Max value (e.g., cloumn 3) for entire CS+ presentation
title('CS+','FontSize',24)
xlabel('Trials')
ylabel('Pupil size')
format_fig;

figure (12131);
plot(1:30,grpstats(all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),3),all_CSminus_avValues_Ext(all_CSminus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Max value (e.g., cloumn 3) for entire CS- presentation
hold on
plot(1:30,grpstats(all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),3),all_CSplus_avValues_Ext(all_CSplus_avValues_Ext(:,1)==uniqueSubIDs_Ext(k),2)))% Max value (e.g., cloumn 3) for entire CS+ presentation
legend('CS-','CS+')
xlabel('Trials')
ylabel('Pupil size')
format_fig;


%TRIAL BY TRIAL FOR OVERALL SCORES FOR MEDIAN, MAX AND MEAN for split of
%image alone and with liquiD:
%MEAN 
figure(1215);
plot(1:30,grpstats(all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),5),all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),2)))
hold on
plot(1:30,grpstats(all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),8),all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),2))) 
legend ('CS-alone','CS-image+liquid')

figure (1216);
plot(1:30,grpstats(all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),5),all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),2))) 
hold on
plot(1:30,grpstats(all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),8),all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),2)))  
legend ('CS+alone','CS+image+liquid')

%MEDIAN for split of image alone and with liquid
figure(1217);
plot(1:30,grpstats(all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),4),all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),2))) 
hold on
plot(1:30,grpstats(all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),7),all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),2))) 

figure (1218);
plot(1:30,grpstats(all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),4),all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),2)))
hold on
plot(1:30,grpstats(all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),7),all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),2)))  

%MAX for split of image alone and with liquid
figure(1219);
plot(1:30,grpstats(all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),3),all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),2))) 
hold on
plot(1:30,grpstats(all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),6),all_CSminus_avValuesE(all_CSminus_avValuesE(:,1)==uniqueSubIDsE(k),2))) 
legend ('CS-alone','CS-image+liquid')

figure (1220);
plot(1:30,grpstats(all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),3),all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),2))) 
hold on
plot(1:30,grpstats(all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),6),all_CSplus_avValuesE(all_CSplus_avValuesE(:,1)==uniqueSubIDsE(k),2))) 
legend ('CS+alone','CS+image+liquid')

%%
%%% SECTION 13: Comparing first and last 5 trials 
figure(1301); 
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
ylim([-100 350])

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
ylim([-100 350])
%%
%%% SECTION 14: Average of each section of trials

figure (1401); 
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

figure (1402);
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

figure (1403); 
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

figure (1404);
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

figure (1405);
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

figure(1406);
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


