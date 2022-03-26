% main function for epoching
% analysis focuses on eye-closed condition
% sampling rate (downsampled): 512Hz
% strategy:
% 1) down-sampling
% 2) epoching
% 3) rej epoch with extrime value 100 uV
% 4) cal power spectrum
% Written by Yamin Li, Feb, 2022 at SJTU




clc;clear;

load('/Volumes/Backup Plus/疼痛数据/pname.mat');

fs = 512;
dt = 1/fs;


pindex = 1;
for p = 1:37
    if p < 20
        pname = ttname{p,1};
        if ismember(pname,{'CDX' 'TYL' 'ZQF' 'CZZ'})
            continue
        end
        filepath = ['/Volumes/Backup Plus/疼痛数据/preprocessed/' pname '_tt.set'];
        pindex = pindex + 1;
        savename = [pname '_tt.set'];
    else
        pname = dzname{p-19,1};
        savename = [pname '_dz.set'];
        pindex = pindex + 1;
        filepath = ['/Volumes/Backup Plus/疼痛数据/preprocessed/' pname '_dz.set'];
    end
    
    

    EEG = pop_loadset(filepath);
    EEG = pop_resample( EEG, fs);
    
    %% epoching
    EEG = pop_epoch( EEG, {  EEG.event.type  }, [0  300]);
    EEG = pop_rmbase( EEG, [],[]);
    epoch_t = 2;
    overlap_t = 1;
    Npepoch = floor(epoch_t/dt);
    Nol = floor(overlap_t/dt);
    EEGdata = EEG.data;  % 64*timepoints
    Ndata = size(EEGdata,2);
    tstart = 1:Nol:(Ndata-Npepoch+1);
    e_index = 1;
    for i = tstart
        data_e(:,:,e_index) = EEGdata(:,i:i+Npepoch-1);
        e_index = e_index + 1;
    end
    
    % reject epoch over 100 uV
    data_e_temp = permute(data_e,[3 2 1]);
    [tr_rej,~] = find((data_e_temp>100) | (data_e_temp < -100));
    tr_rej = unique(tr_rej);
    data_e(:,:,tr_rej) = [];
    
    if p < 20
        groupindex = 1;
        EpochData{p,groupindex} = data_e;
    else
        groupindex = 2;
        EpochData{p-19,groupindex} = data_e;
    end
    clear data_e data_e_temp tr_rej EEGdata
    
end