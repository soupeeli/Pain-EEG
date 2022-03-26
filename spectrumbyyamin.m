% Codes for calculating power spectrum and band power
%
% Copyright(c) 2020 Yamin Li, School of Biomedical Engineering
% Shanghai Jiao Tong University

clc
clear all;
close all;

%% loading data and preparation

% 导入dataset
EEG=pop_loadset('E:\liyamin\2020615renliyuan eeg\social2 better\20200615social_epoch.set');


Ne=size(EEG.data,1); % 导联数
Np=size(EEG.data,2); % 每个epoch点数
Nt=size(EEG.data,3); % trial数
EEGdata=EEG.data;

% define ROI: PL-posterior left; PR-posterior right; PM-posterior middle
chan_id_F = []; 
chan_id_LT = [];
chan_id_RT = []; 
chan_id_C = [];
chan_id_PO =[];

ROI_F = {'Fp1','Fp2','Fpz','AF3','AF7','AF4','AF8'...
    'F7','F5','F3','F1','Fz','F2','F4','F6','F8'};
for i = 1:length(ROI_F)
    chan_id_F = [chan_id_F find(strcmpi(ROI_F{i},{EEG.chanlocs.labels}))];
end

ROI_LT = {'FT7','FC5','FC3', ...
    'T7','C5','C3', ...
    'TP7','CP5','CP3'};
for i = 1:length(ROI_LT)
    chan_id_LT = [chan_id_LT find(strcmpi(ROI_LT{i},{EEG.chanlocs.labels}))];
end

ROI_RT = {'FT8','FC6','FC4', ...
    'T8','C6','C4', ...
    'TP8','CP6','CP4'};
for i = 1:length(ROI_RT)
    chan_id_RT = [chan_id_RT find(strcmpi(ROI_RT{i},{EEG.chanlocs.labels}))];
end

ROI_C = {'FC1','FCz','FC2',...
    'C1','Cz','C2', ...
    'CP1','CP2'};
for i = 1:length(ROI_C)
    chan_id_C = [chan_id_C find(strcmpi(ROI_C{i},{EEG.chanlocs.labels}))];
end

ROI_PO = {'P3','P5','P7','PO3','PO7','O1', ...
    'P4','P6','P8','PO4','PO8','O2', ...
    'P1','P2','Pz','Oz','POz','O1','O2'};
for i = 1:length(ROI_PO)
    chan_id_PO = [chan_id_PO find(strcmpi(ROI_PO{i},{EEG.chanlocs.labels}))];
end


ques='Please input ROI: all/F/LT/RT/C/PO standing for all channels/frontal/left/right temporal/central/parietal occipital, respectively: \n';
ROI=input(ques);

%% FFT
% EEG data format: electrodes * points * trials

fs=1000;
N=Np;
n=0:N-1;
nfft=2^(nextpow2(N));
f_fft=(0:nfft/2-1)/nfft*fs;

% plot with grid
h1=figure;                                       

for i = 1:Nt
    data=EEGdata(:,:,i);
    mag=abs(fft(data,nfft,2)); % fast fourier transform
    Pxx(:,:,i)=mag.^2/N/fs*2;        
%         plot(20*log(power_spectrum)); hold on; grid on;   %frequency
end
Pxx_AVG=mean(Pxx,3);

%% pwelch method
win=2000;
overlap=1800;


for i = 1:Nt
    for channel = 1:57
        [Pwelch(channel,:,i),f_welch] = pwelch(EEGdata(channel,:,i),win,overlap,[],fs);
    end
end
    
% 数据段平均
Pwelch_AVG=mean(Pwelch,3);



%% plot
switch ROI
    case 'all'
        Pxx_AVG_ROI=mean(Pxx_AVG,1);
        Pwelch_AVG_ROI=mean(Pwelch_AVG,1);
    case 'F'
        Pxx_AVG_ROI=mean(Pxx_AVG(chan_id_F,:),1);
        Pwelch_AVG_ROI=mean(Pwelch_AVG(chan_id_F,:),1);
    case 'LT'
        Pxx_AVG_ROI=mean(Pxx_AVG(chan_id_LT,:),1);
        Pwelch_AVG_ROI=mean(Pwelch_AVG(chan_id_LT,:),1);
    case 'RT'
        Pxx_AVG_ROI=mean(Pxx_AVG(chan_id_RT,:),1);
        Pwelch_AVG_ROI=mean(Pwelch_AVG(chan_id_RT,:),1);
    case 'C'
        Pxx_AVG_ROI=mean(Pxx_AVG(chan_id_C,:),1);
        Pwelch_AVG_ROI=mean(Pwelch_AVG(chan_id_C,:),1);
    case 'PO'
        Pxx_AVG_ROI=mean(Pxx_AVG(chan_id_PO,:),1);
        Pwelch_AVG_ROI=mean(Pwelch_AVG(chan_id_PO,:),1);
end


y_fft =  Pxx_AVG_ROI(1:nfft/2);
y_welch = Pwelch_AVG_ROI;

%     plot(f,10*log(y{s})); hold on; grid on;
plot(f_fft,y_fft); hold on; grid on;
plot(f_welch,y_welch,'blue');

set(gca, 'xlim', [2 40], 'fontsize',15);%, 'ylim', [0 3]
xlabel('Frequence(Hz)','FontSize', 15);
ylabel('Amplitude(μV^2)','FontSize', 15);
title('Power Spectrum','FontSize', 15);

fprintf('结果存在变量y里，横坐标频率在变量f里\n');

% close(h1)

%% calculate band power
delta=find(f_welch>0.1 & f_welch < 4);
theta=find(f_welch>4 & f_welch < 8);
alpha=find(f_welch>8 & f_welch < 12);
beta=find(f_welch>12 & f_welch < 30);
gamma=find(f_welch>30 & f_welch < 40);

power_delta=sum(y_welch(delta));
power_theta=sum(y_welch(theta));
power_alpha=sum(y_welch(alpha));
power_beta=sum(y_welch(beta));
power_gamma=sum(y_welch(gamma));
fprintf('delta power=%f\n',power_delta);
fprintf('theta power=%f\n',power_theta);
fprintf('alpha power=%f\n',power_alpha);
fprintf('beta power=%f\n',power_beta);
fprintf('gamma power=%f\n',power_gamma);
