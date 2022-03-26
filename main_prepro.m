% main function for preprocessing the pain data 
% Feb. 2022
% - Yamin
clc;clear;

load('/Volumes/Backup Plus/疼痛数据/pname.mat');
for p = 22:22
    if p < 20
        pname = ttname{p,1};
        event = ttevent{p,1};
        if ismember(pname,{'CDX' 'TYL' 'ZQF'})
            continue
        end
        savename = [pname '_tt.set'];
        filepath = ['/Volumes/Backup Plus/疼痛数据/mxtt/' pname '.edf'];
    else
        pname = dzname{p-19,1};
        event = dzevent{p-19,1};
        savename = [pname '_dz.set'];
        filepath = ['/Volumes/Backup Plus/疼痛数据/jkdz/' pname '.edf'];
    end
    
    EEG = pop_biosig(filepath);
    EEG = pop_select( EEG, 'nochannel',{'EMG1','EMG2','EOG EOGL','EOGR','MK'});
    EEG = pop_eegfiltnew(EEG, 'locutoff',1);
    EEG = pop_eegfiltnew(EEG, 'hicutoff',100);
    EEG = pop_eegfiltnew(EEG, 'locutoff',45,'hicutoff',55,'revfilt',1);
    EEG = pop_epoch( EEG, {event}, [-304  304], 'newname', 'EDF file epochs', 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, [],[]);
    
    EEG = pop_saveset( EEG, 'filename',savename,'filepath','/Volumes/Backup Plus/疼痛数据/preprocessed/');
end


