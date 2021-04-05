%% Model Script to analyse Anti-Saccade Eye-Link data
% author: thomas.andrillon@monash.edu
clear all
close all

run localdef.m

addpath(genpath(lscpTools_path));
addpath(edf2mat_path);

files=dir([eyet_path filesep '*' filesep '*.edf']);

%% Loop on files
redo=1;
for n=1:length(files)
    subID=files(n).name(1:end-4);
    savename=['Kavaliotis_' subID];
    
    if exist([files(n).folder filesep savename '.mat'])==0 || redo==1
        fprintf('... converting %s from EDF to .mat file\n',subID)
        % import into matlab
        myedf=Edf2Mat([files(n).folder filesep subID '.edf']);
        
        % clean from useless info
        %     data=myedf;
        %     list_fields_toremove={'AUTHOR','AUTHOREMAIL','COPYRIGHT','VERSION','VERSIONDATE','CHANGELOG',...
        %         'RECORDING_STATES',};
        %     data=rmfiled(data,'AUTHOR');
        EL_headers=[];
        EL_headers=myedf.Header;
        EL_headers.Fs=unique(diff(myedf.timeline))*1000'; % in Hertz
        
        EL_data=[];
        EL_data.time=myedf.Samples.time;
        EL_data.pupilSize=myedf.Samples.pupilSize;
        EL_data.posX=myedf.Samples.posX;
        EL_data.posY=myedf.Samples.posY;
        
        EL_events=[];
        EL_events.Events.time=myedf.Events.Messages.time;
        EL_events.Events.type=myedf.Events.Messages.info;
        EL_events.StartRec=myedf.Events.Start.time;
        EL_events.EndRec=myedf.Events.End.time;
        
        EL_events.Fix.start=myedf.Events.Efix.start;
        EL_events.Fix.end=myedf.Events.Efix.end;
        EL_events.Fix.duration=myedf.Events.Efix.duration;
        EL_events.Fix.posX=myedf.Events.Efix.posX;
        EL_events.Fix.posY=myedf.Events.Efix.posY;
        EL_events.Fix.pupilSize=myedf.Events.Efix.pupilSize;
        
        EL_events.Blinks.start=myedf.Events.Eblink.start;
        EL_events.Blinks.end=myedf.Events.Eblink.end;
        EL_events.Blinks.duration=myedf.Events.Eblink.duration;
        
        
        EL_events.Sacc.start=myedf.Events.Esacc.start;
        EL_events.Sacc.end=myedf.Events.Esacc.end;
        EL_events.Sacc.duration=myedf.Events.Esacc.duration;
        EL_events.Sacc.posX_start=myedf.Events.Esacc.posX;
        EL_events.Sacc.posY_start=myedf.Events.Esacc.posY;
        EL_events.Sacc.posX_end=myedf.Events.Esacc.posXend;
        EL_events.Sacc.posY_end=myedf.Events.Esacc.posYend;
        EL_events.Sacc.velo=myedf.Events.Esacc.pvel;
        
        save([files(n).folder filesep savename],'EL_headers','EL_data','EL_events');
    else
        fprintf('... load %s from .mat file\n',subID)
        load([files(n).folder filesep savename])
    end
       
    %     ). Blinks were linearly interpolated from 200 ms before to 200 ms
    % 756 after automatically identified blinks, and the interpolated pupil data was then low-pass filtered (< 6
    % 757 Hz, second order butterworth).
    [data_pupil,  filt_pupilSize]=get_EyeLink_cleanpupil(EL_data.pupilSize,EL_headers.Fs,EL_data.time,EL_events);
    EL_data.clean_pupilSize=data_pupil;
    EL_data.filt_pupilSize=filt_pupilSize;
    save([files(n).folder filesep savename '_clean'],'EL_headers','EL_data','EL_events');
    %     else
    %         fprintf('... %s already imported\n',subID)
    % end
end

%%
figure;
plot(EL_data.pupilSize);
hold on;
plot(EL_data.clean_pupilSize);
plot(EL_data.filt_pupilSize);
legend({'raw','blink-corrected','filtered'})
format_fig;
xlabel('sample')
ylabel('pupil size')
