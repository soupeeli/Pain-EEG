


clc;clear;

load('/Volumes/Backup Plus/疼痛数据/pname.mat');

fs = 512;

%% load the required geometrical information
        load('standard_bem.mat');
        load('elec_aligned.mat');
        headmodel_eeg = vol;
        
        headmodel_eeg = ft_convert_units(headmodel_eeg, 'mm');

filepath = '/Volumes/Backup Plus/疼痛数据/source trial/';
savepath = '/Volumes/Backup Plus/疼痛数据/conn result/';

%% Start calculation
for p = 1:19
    if p < 20
        pname = ttname{p,1};
        if ismember(pname,{'CDX' 'TYL' 'ZQF' 'CZZ'})
            continue
        end
        filename = [pname '_tt_oritr'];
        savename = [pname '_tt_conn_'];
    else
        pname = dzname{p-19,1};
        filename = [pname '_dz_oritr'];
        savename = [pname '_dz_conn_'];
    end
    load([filepath filename '.mat']);  % parameter: source
    
    for fre = 1:5
        %% Setting the final filepath
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
    
        % filtering the data
        % frequency band of interest:δ波（1-4Hz），θ波（4-8Hz），α波（8-13Hz），β波（13-30Hz），γ波（30-100Hz，对45~55Hz进行带通滤波后）
        cfg = [];
        cfg.hpfilter       = 'yes';        % enable high-pass filtering
        cfg.lpfilter       = 'yes';        % enable low-pass filtering
        cfg.hpfreq         = flow;           % set up the frequency for high-pass filter
        cfg.lpfreq         = fhigh;          % set up the frequency for low-pass filter
        data_filtered      = ft_preprocessing(cfg,data_clean);
        
        Ntr = size(data_filtered.trial,2);
        data_tr= [];
        for tr = 1:Ntr
            data_tr(:,:,tr) = data_filtered.trial{1,tr};
        end
        
        
        %% extracting instaneous phase
        data_complex = [];
        for ch = 1:size(data_tr,1)
            fprintf('calculating phase %d/%d\n', ch, size(data_tr,1));
            data_complex(ch,:,:) = hilbert(squeeze(data_tr(ch, :, :)));
%             Ang_source(ch, :, :) = angle(source_complex(ch,:,:));
        end
        %% calculating phase-based brain network
        
        % --- PLV - phase locking value ---
        [PLV.fullconn,PLI.fullconn,wPLI.fullconn] = ...
            calphasesyn(data_complex);
        
        %% Thresholding & calculating the graph index
        % Thresholding
        thres = 0.1;
        PLV.conn_thre = threshold_proportional(PLV.fullconn,thres);
        PLV.conn_thre_b = PLV.conn_thre;
        PLV.conn_thre_b(PLV.conn_thre_b ~=0)=1;
        
        PLI.conn_thre = threshold_proportional(PLI.fullconn,thres);
        PLI.conn_thre_b = PLI.conn_thre;
        PLI.conn_thre_b(PLI.conn_thre_b ~=0)=1;
        
        wPLI.conn_thre = threshold_proportional(wPLI.fullconn,thres);
        wPLI.conn_thre_b = wPLI.conn_thre;
        wPLI.conn_thre_b(wPLI.conn_thre_b ~=0)=1;
        
        % Local graph index: degree, local efficiency
        
        % degree
        PLV.degree = degrees_und(PLV.conn_thre_b);
        PLI.degree = degrees_und(PLI.conn_thre_b);
        wPLI.degree = degrees_und(wPLI.conn_thre_b);
        
        % local efficiency
        PLV.loceff = efficiency_bin(PLV.conn_thre_b,1);
        PLI.loceff = efficiency_bin(PLI.conn_thre_b,1);
        wPLI.loceff = efficiency_bin(wPLI.conn_thre_b,1);
        
        
        % 4 Global graph index: global efficiency; small-worldness; 
        %                       global clustering coefficient; char path;
        
        % global efficiency
        PLV.glbeff = efficiency_bin(PLV.conn_thre_b);
        PLI.glbeff = efficiency_bin(PLI.conn_thre_b);
        wPLI.glbeff = efficiency_bin(wPLI.conn_thre_b);
        
        % charpath
        plvD_temp = distance_bin(PLV.conn_thre_b);
        PLV.chpath = charpath(plvD_temp,0,0);
        
        pliD_temp = distance_bin(PLI.conn_thre_b);
        PLI.chpath = charpath(pliD_temp,0,0);
        
        wpliD_temp = distance_bin(wPLI.conn_thre_b);
        wPLI.chpath = charpath(wpliD_temp,0,0);
        
        % small worldness
        PLV.smwolrd = calsmallworld(PLV.conn_thre_b,1);
        PLI.smwolrd = calsmallworld(PLI.conn_thre_b,1);
        wPLI.smwolrd = calsmallworld(wPLI.conn_thre_b,1);
        

        % global CC
        PLV.glbCC = mean(clustering_coef_bu(PLV.conn_thre_b));
        PLI.glbCC = mean(clustering_coef_bu(PLI.conn_thre_b));
        wPLI.glbCC = mean(clustering_coef_bu(wPLI.conn_thre_b));
        
        
        
        save([savepath savename frename '.mat'],'PLV','PLI','wPLI');
        
        clear data_complex  PLV PLI wPLI
        fprintf(['person ' num2str(p) 'fre ' num2str(fre) ' finished!\n'])
    end
    clear data_clean

end

        
               


%% Functions


function [PLV, PLI, wPLI] = calphasesyn(complexdata)

Nchan = size(complexdata, 1);
Ntr = size(complexdata, 3);

PLV = zeros(Nchan,Nchan);
PLI = zeros(Nchan,Nchan);
wPLI = zeros(Nchan,Nchan);

Ang = angle(complexdata);

for i = 1:Nchan
    for j = i+1:Nchan
        % PLV
        deltaphase = squeeze(Ang(i, :, :)) - squeeze(Ang(j, :, :));
        PLV(i, j) = sum(abs(mean(exp(1i*deltaphase), 1))) / Ntr; % checked: no problem 2021/10/13
        
        % Initial PLI, wPLI
        PLI_temp = 0;
        wPLI_temp = 0;
        for tr = 1:Ntr
            channelData = squeeze(Ang(i, :, tr));
            compareChannelData = squeeze(Ang(j, :, tr));
            phasediff = channelData - compareChannelData;
            
            % PLI
            v=sign(sin(phasediff));
            PLI_temp = PLI_temp + abs(mean(v));
            
            % wPLI
            img = angle(complexdata(i, :, tr)./complexdata(j, :, tr));
            wPLI_temp = wPLI_temp + abs(sum(abs(img).*sign(img))/sum(abs(img)));            
        end
        PLI(i,j) = PLI_temp/Ntr;
        wPLI(i,j) = wPLI_temp/Ntr;
        
        PLI(j,i) = PLI(i,j);
        wPLI(j,i) = wPLI(i, j);
        PLV(j, i) = PLV(i, j);
    end
    fprintf(['channel ' num2str(i) ' finished!\n']);
end

end


%%
function [plv] = eegPLV_phase(Data)
% Computes the Phase Locking Value (PLV) for an EEG dataset.
% 
%
% Input parameters:
%   eegData is a 3D matrix numChannels x numTimePoints x numTrials
%   srate is the sampling rate of the EEG data
%   toi is the time of interest for analysis. form e.g.: [300 100] for
%   300-1000 points of the eeg sequences.
%
% Output parameters:
%   plv is a 2D matrix - 
%     numChannels x numChannels
%--------------------------------------------------------------------------
% Example: Consider a 28 channel EEG data sampled @ 1000 Hz with 231 trials,
% where each trial lasts for 1 seconds. Below is an example of how to
% use this function.
%
%   eegData = rand(28, 1000, 231);
%   toi = [500 1000]; % only calculating the 500-1000 points
%   [plv] = eegPLV(eegData, toi);
% 
%   NOTE:
% 
% Also note that in order to extract the PLV between channels 17 and 20, 
% use plv(:, 17, 20, :) and NOT plv(:, 20, 17, :). The smaller channel 
% number is to be used first.
%--------------------------------------------------------------------------
% 
% Reference:
%   Lachaux, J P, E Rodriguez, J Martinerie, and F J Varela. 
%   Measuring phase synchrony in brain signals.? 
%   Human brain mapping 8, no. 4 (January 1999): 194-208. 
%   http://www.ncbi.nlm.nih.gov/pubmed/10619414.
% 
%--------------------------------------------------------------------------
% Written by: 
% Yamin Li
% Neuroengineering Lab, SJTU
% 01 Dec 2020



Nchan = size(Data, 1);
numTrials = size(Data, 3);
Ang = Data;

plv = zeros(Nchan,Nchan);

for i = 1:Nchan
    channelData = squeeze(Ang(i, :, :));
    for j = i+1:Nchan
        compareChannelData = squeeze(Ang(j, :, :));
        phasediff = channelData - compareChannelData;
        plv(i, j) = sum(abs(mean(exp(1i*phasediff), 1))) /numTrials; % checked: no problem 2021/10/13
        plv(j, i) = plv(i, j);
    end
end
end



%%
function [pli] = PLI_phase(phaseData)
% channel * time * trial


Nchan = size(phaseData, 1);
numTrials = size(phaseData, 3);
Ang = phaseData;
pli = zeros(Nchan,Nchan);

for i = 1:Nchan
    channelData = squeeze(Ang(i, :, :));
    for j = i+1:Nchan
        compareChannelData = squeeze(Ang(j, :, :));
        phasediff = channelData - compareChannelData;
        sum_sign = [];
        for tr = 1:numTrials
            delta_p = phasediff(:,tr);
            v=sign(sin(delta_p));
            sum_sign=sum_sign+abs(mean(v));
        end
        pli(i,j) = sum_sign/numTrials;
        pli(j,i) = pli(i,j);
    end
end


end

%% wPLI
function [wPLI] = cal_wPLI(complexdata)

Nchan = size(complexdata, 1);
Ntr = size(complexdata, 3);
wPLI = zeros(Nchan,Nchan);

for i = 1:Nchan
    for j = i+1:Nchan
        wPLI_temp = [];
        for tr = 1:numTrials
            img = angle(complexdata(i, :, tr)./complexdata(j, :, tr));
            wPLI_temp(i, j,tr) = abs(sum(abs(img).*sign(img))/sum(abs(img)));            
        end
        wPLI(i,j) = mean(wPLI_temp,3);
        wPLI(j,i) = wPLI(i, j);
    end
end


end

%%
function S = calsmallworld(A,flag)

% get its basic properties
n = size(A,1);  % number of nodes
k = sum(A);  % degree distribution of undirected network
m = sum(k)/2;
K = mean(k); % mean degree of network

%% computing small-world-ness using the analytical approximations for the E-R graph

[expectedC,expectedL] = ER_Expected_L_C(K,n);  % L_rand and C_rand

[S,~,~] = small_world_ness(A,expectedL,expectedC,flag);  % Using WS clustering coefficient

end

