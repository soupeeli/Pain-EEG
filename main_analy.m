% Script for preprocessing and epoching pain data using fieldtrip toolbox

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
savepath = '/Volumes/Backup Plus/疼痛数据/source trial/';

fs = 512;
pindex = 1;

%% load the required geometrical information
        load('standard_bem.mat');
        load('elec_aligned.mat');
%         load('standard_sourcemodel3d10mm.mat');
        headmodel_eeg = vol;
        
        headmodel_eeg = ft_convert_units(headmodel_eeg, 'mm');
%         sourcemodel = ft_convert_units(sourcemodel,'mm');


%% Start calculation
for p = 8:8
    if p < 20
        pname = ttname{p,1};
        event = ttevent{p,1};
        if ismember(pname,{'CDX' 'TYL' 'ZQF' 'CZZ'})
            continue
        end
        filepath = ['/Volumes/Backup Plus/疼痛数据/preprocessed chan/' pname '_tt.set'];
        pindex = pindex + 1;
        savename_s = [pname '_tt_s_'];
        savename_oritr = [pname '_tt_oritr'];
    else
        pname = dzname{p-19,1};
        event = dzevent{p-19,1};
        savename_s = [pname '_dz_s_'];
        savename_oritr = [pname '_dz_oritr'];
        pindex = pindex + 1;
        filepath = ['/Volumes/Backup Plus/疼痛数据/preprocessed chan/' pname '_dz.set'];
    end
    
    cfg = []; 
    cfg.dataset = filepath; 
    cfg.channel = 'all'; % this is the default
    cfg.reref = 'yes';
    cfg.refmethod = 'avg';
    cfg.refchannel = 'all';
    data = ft_preprocessing(cfg);
    
    %% epoching and rejecting artifacts
    cfg         = [];
    cfg.length  = 2;
    cfg.overlap = 0.5;
    data        = ft_redefinetrial(cfg, data);

    % this step is needed to 1) remove the DC-component, and to 2) get rid of a few segments of data at
    % the end of the recording, which contains only 0's.
    cfg        = [];
%     cfg.demean = 'yes';
    cfg.trials = 1:(numel(data.trial)-6);
    data       = ft_preprocessing(cfg, data);
    
    % reject epochs with artifacts
    cfg = [];
    cfg.continuous = 'no';
    cfg.artfctdef.threshold.min       = -100;
    cfg.artfctdef.threshold.max       = 100;
    [cfg, artifact] = ft_artifact_threshold(cfg, data);
    data_clean = ft_rejectartifact(cfg, data);
    clear data
    
%% visualize the coregistration of sensors, headmodel, and sourcemodel.
    
% % %     figure(1);
% % %     % make the headmodel surface transparent
% % % 
% % %     ft_plot_headmodel(headmodel_eeg, 'edgecolor', 'none'); 
% % %     ft_plot_sens(data_filtered.elec);
% % %     view([45 -15 0])
% % % 
% % %     cfg = []; 
% % %     cfg.method = 'interactive'; 
% % %     cfg.headshape = headmodel_eeg.bnd(1); 
% % %     cfg.elec = data_filtered.elec; 
% % %     elec_aligned = ft_electroderealign(cfg);

    data_clean.elec = elec_aligned; % make sure the aligned electrodes are updated

% % %     % examine the location
% % %     ft_plot_headmodel(headmodel_eeg, 'edgecolor', 'none'); 
% % %     ft_plot_sens(data_clean.elec);
% % %     view([45 -15 0]);

    %% prepare source model
    % Compute Grid
    cfg                 = [];
    cfg.elec            = elec_aligned;
    cfg.headmodel       = headmodel_eeg;
    cfg.resolution = 10;
    cfg.unit       = 'mm';
    cfg.inwardshift     = -1.5;
    sourcemodel         = ft_prepare_sourcemodel(cfg);

    %% prepare lead field
    cfg         = [];
    cfg.elec    = elec_aligned;            
    cfg.channel = data_clean.label;  
    cfg.sourcemodel    = sourcemodel;           % 2002v source points
    cfg.headmodel = headmodel_eeg;                      % volume conduction model
    leadfield     = ft_prepare_leadfield(cfg);
        
    
    
    % *** important ***
    % If you would know that the subsequent analysis would be limited to a specific frequency range in the data (e.g., everything above 30 Hz), you could first apply a filter using ft_preprocessing (e.g., cfg.hpfilter=yes and cfg.hpfreq=30) prior to computing the covariance and the spatial filter.
    % *** --------- ***
    
    % filtering the data
    % frequency band of interest:δ波（1-4Hz），θ波（4-8Hz），α波（8-13Hz），β波（13-30Hz），γ波（30-100Hz，对45~55Hz进行带通滤波后）
    for fre = 1:5
        switch fre
            case 1
                flow = 1;
                fhigh = 4;
            case 2
                flow = 4;
                fhigh = 8;
            case 3
                flow = 8;
                fhigh = 13;
            case 4
                flow = 13;
                fhigh = 30;
            case 5
                flow = 30;
                fhigh = 100;
        end
        cfg = [];
        cfg.hpfilter       = 'yes';        % enable high-pass filtering
        cfg.lpfilter       = 'yes';        % enable low-pass filtering
        cfg.hpfreq         = flow;           % set up the frequency for high-pass filter
        cfg.lpfreq         = fhigh;          % set up the frequency for low-pass filter
        data_filtered      = ft_preprocessing(cfg,data_clean);
        
        
        %% calculate covariance
        cfg                   = [];
        cfg.covariance        = 'yes';
        cfg.preproc.demean = 'yes';
        cfg.covariancewindow  = 'all';
        cfg.keeptrials        = 'yes';
        tlock = ft_timelockanalysis(cfg, data_filtered);
        
        
        %% source analysis
        
        cfg                  = [];
        cfg.method           = 'lcmv';
        cfg.sourcemodel      = leadfield; % leadfield
        cfg.headmodel        = headmodel_eeg; % volume conduction model (headmodel)
        cfg.lcmv.keepfilter  = 'yes';
        cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
        cfg.lcmv.lambda       = '5%';  %?regularization parameter of 5%
%         cfg.keeptrials        = 'yes';
%         cfg.rawtrial = 'yes';
        source_tmp               = ft_sourceanalysis(cfg, tlock); 
        filters               = cell2mat(source_tmp.avg.filter(source_tmp.inside)); % keep filters

        %% combine filters with raw data to "beam" scalp data to source space
        source = [];
        % for each trial, apply filters to the recorded data
        for j = 1:size(data_filtered.trial,2)
            source.trial(:,:,j) = single(filters*squeeze(data_filtered.trial{1,j}));
        end
        
        switch fre
            case 1
                frename = 'delta';
            case 2
                frename = 'theta';
            case 3
                frename = 'alpha';
            case 4
                frename = 'beta';
            case 5
                frename = 'gamma';
        end
        
        savename_final = [savename_s frename '.mat'];
        save([savepath savename_final],'source');
        fprintf(['person ' num2str(p) ' frequency ' num2str(fre) ' finished!\n']);
        clear data_filtered
    end
    save([savepath savename_oritr '.mat'],'data_clean')
    clear source data_clean 

    % move the DC component of the clean trial. if this step is placed
    % before rejecting the artifact, there would be an error
%     cfg = [];
%     cfg.demean = 'yes';
%     data_clean       = ft_preprocessing(cfg, data_clean);
    
    
    %The LCMV spatial filter is computed here without applying any time-domain 
    %filters. Consequently, it will have to suppress all noise in the data in all frequency bands. 
    %The spatial filter derived from the broadband data allows us to compute a broadband source level time series.

     
    
end

sourceinfo            = [];
sourceinfo.sampleinfo = data_filtered.sampleinfo;
sourceinfo.time       = data_filtered.time;
sourceinfo.fsample    = 512;
sourceinfo.pos = source_tmp.pos(source_tmp.inside,:);
% create labels for each virtual electrode
for c = 1 : sum(sourcemodel.inside)
    label{c,1} = ['S' num2str(c)];
end
sourceinfo.label = label;
save([savepath 'sourceinfo.mat'],'sourceinfo');
