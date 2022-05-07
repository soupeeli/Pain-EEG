% temp script to edit the channel location information
clc;clear;


load('/Volumes/Backup Plus/疼痛数据/pname.mat');

fs = 512;

pindex = 1;
for p = 1:37
    if p < 20
        pname = ttname{p,1};
        event = ttevent{p,1};
        if ismember(pname,{'CDX' 'TYL' 'ZQF' 'CZZ'})
            continue
        end
        filepath = ['/Volumes/Backup Plus/疼痛数据/preprocessed/' pname '_tt.set'];
        pindex = pindex + 1;
        savename = [pname '_tt.set'];
    else
        pname = dzname{p-19,1};
        event = dzevent{p-19,1};
        savename = [pname '_dz.set'];
        pindex = pindex + 1;
        filepath = ['/Volumes/Backup Plus/疼痛数据/preprocessed/' pname '_dz.set'];
    end
    
    EEG = pop_loadset(filepath);
    Nchan = size(EEG.chanlocs,1);
    for c = 1:Nchan
        EEG.chanlocs(c).labels = erase(EEG.chanlocs(c).labels,'EEG ');
    end
    EEG=pop_chanedit(EEG, 'lookup','/Users/soupee/Documents/MATLAB/eeglab2021.1/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
    EEG = pop_epoch( EEG, {  EEG.event.type  }, [0  300]);
    EEG = pop_rmbase( EEG, [],[]);
    EEG = pop_resample( EEG, fs);
    EEG = pop_saveset( EEG, 'filename',savename,'filepath','/Volumes/Backup Plus/疼痛数据/preprocessed chan/');

    
end