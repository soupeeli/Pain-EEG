% main function for extracting individual dominant peak frequency
% Code by Yamin Li, Shanghai Jiao Tong University
clear;clc;close all;

load('/Volumes/Backup Plus/ÌÛÍ´Êý¾Ý/powerresult.mat');
flow = 6;
fhigh = 14;
frange = find(f_welch>=6 & f_welch <=14);
f_seq = f_welch(frange);

%% preparing the data
power_dz_all = squeeze(mean(power_dz,1));
power_tt_all = squeeze(mean(power_tt,1));

%% calculation

pnum_dz = size(power_dz_all,2);
pnum_tt = size(power_tt_all,2);

for p = 1:pnum_tt
    [pks,locs] = findpeaks(power_tt_all(frange,p),'sortstr','descend');
    pks_tt(p) = pks(1);
    locs_tt(p) = locs(1);
    fre_tt(p) = f_seq(locs(1));
end

for p = 1:pnum_dz
    [pks,locs] = findpeaks(power_dz_all(frange,p),'sortstr','descend');
    pks_dz(p) = pks(1);
    locs_dz(p) = locs(1);
    fre_dz(p) = f_seq(locs(1));
end


pks_tt_dB = 10*log10(pks_tt);
pks_dz_dB = 10*log10(pks_dz);

%% Plot
figure

scatter(fre_tt,pks_tt_dB,50,[0.8500 0.3250 0.0980],'+','LineWidth',1.5);
hold on;
scatter(fre_dz,pks_dz_dB,50,[0.4660 0.6740 0.1880],'o','LineWidth',1.5);
legend('Pain','Control');
xlim([6 14]);
ylim([-5 25]);
xlabel('Frequency (Hz)', 'FontSize',12);
ylabel('Power 10*log10 (\muV^2/Hz)', 'FontSize',12);


