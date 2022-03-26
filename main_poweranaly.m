% main function for conducting spectural analysis

% load('/Volumes/Backup Plus/疼痛数据/Epochdata_dz.mat')
% load('/Volumes/Backup Plus/疼痛数据/Epochdata_tt.mat')
clear;clc;
fs = 512;
f_theta = [4 8];
f_alpha = [8 14];
f_beta = [14 30];
f_gamma = [60 99];


% % % for group = 1:2
% % %     if group == 1
% % %         pnum = 19;
% % %     else
% % %         pnum = 18;
% % %     end
% % %     pindex = 1;
% % %     for p = 1:pnum
% % %         if group == 1
% % %             epochdata = EpochData_tt{p,1};
% % %         else
% % %             epochdata = EpochData_dz{p,1};
% % %         end
% % %         if isempty(epochdata)
% % %             continue;
% % %         end
% % %         Nt = size(epochdata,3);
% % %         for i = 1:Nt
% % %             for channel = 1:64
% % %                 [Pwelch(channel,:,i),f_welch] = pwelch(epochdata(channel,:,i),1024,0,[],fs);
% % %             end
% % %         end
% % %         Pwelch_AVG=mean(Pwelch,3);
% % %         if group == 1
% % %             power_tt(:,:,pindex) = Pwelch_AVG;
% % %         else
% % %             power_dz(:,:,pindex) = Pwelch_AVG;
% % %         end
% % %         pindex = pindex + 1;
% % %         clear Pwelch Pwelch_AVG
% % %         fprintf(['Group ' num2str(group) ', p ' num2str(p) ' finished\n']);
% % %     end
% % % end
% % %     
% 数据段平均


%% Topographical t test

ifnorm = 0;


load('/Volumes/Backup Plus/疼痛数据/powerresult.mat')
EEG = pop_loadset('/Volumes/Backup Plus/疼痛数据/64chanloc.set');
chanloc = EEG.chanlocs;


% find band point
frange{1} = find(f_welch>=f_theta(1) & f_welch<f_theta(2));
frange{2} = find(f_welch>=f_alpha(1) & f_welch<f_alpha(2));
frange{3} = find(f_welch>=f_beta(1) & f_welch<f_beta(2));
frange{4} = find(f_welch>=f_gamma(1) & f_welch<f_gamma(2));
frange_all = find(f_welch>=1 &f_welch<100);

figure
for fre = 1:4
    frange_temp = frange{fre};
    power_tt_temp = sum(power_tt(:,frange_temp,:),2);
    power_dz_temp = sum(power_dz(:,frange_temp,:),2);
    
    
    power_tt_fulltemp = sum(power_tt(:,frange_all,:),2);
    power_tt_rel = power_tt_temp./power_tt_fulltemp;
    power_dz_fulltemp = sum(power_dz(:,frange_all,:),2);
    power_dz_rel = power_dz_temp./power_dz_fulltemp;
    
    % spectra
    for f = 1:513
        spe_tt_rel(:,f,:) = power_tt(:,f,:)./power_tt_fulltemp;
        spe_dz_rel(:,f,:) = power_dz(:,f,:)./power_dz_fulltemp;
    end
    
    
    for chan = 1:64
        if ifnorm == 1
            power_tt_chantemp = power_tt_rel(chan,:);
            power_dz_chantemp = power_dz_rel(chan,:);
        else
            power_tt_chantemp = power_tt_temp(chan,:);
            power_dz_chantemp = power_dz_temp(chan,:);
        end
%         power_tt_chantemp = power_tt_temp(chan,:)./power_tt_fulltemp(chan,:);
%         power_dz_chantemp = power_dz_temp(chan,:)./power_dz_fulltemp(chan,:);

        [h{fre}(chan,:),p{fre}(chan,:),~,stats] = ttest2(power_tt_chantemp,power_dz_chantemp);
        T{fre}(chan,:) = stats.tstat;
    end
    
    
    if fre == 2 && ifnorm == 1
        cup = 0.4;
    elseif ifnorm == 1
        cup = 0.2;
    elseif fre == 2
        cup = 110;
    else
        cup = 80;
    end
    
% %     subplot(4,4,fre)
% %     topoplot(T{fre},chanloc,'style','map');
% %     caxis([-4 4])
% % %     colorbar
% %     
% %     subplot(4,4,fre+4)
% %     topoplot(p{fre},chanloc,'style','map','conv','off');
% %     colormap(flipud(jet))
% %     caxis([0 0.05]);

    
    subplot(2,4,fre)
    topoplot(mean(power_tt_temp,3),chanloc,'style','map');
    caxis([0 cup])
    colorbar
    
    subplot(2,4,fre+4)
    topoplot(mean(power_dz_temp,3),chanloc,'style','map');
    caxis([0 cup])
    colorbar

    
    
end

% plot p value
figure
for fre = 1:4
    subplot(1,4,fre)
    topoplot(p{fre},chanloc,'style','map','conv','off');
    colormap(flipud(parula))
    caxis([0 0.05]);
end

% plot power value
figure
for fre = 1:4
    subplot(2,4,fre)
    topoplot(mean(power_tt_rel,3),chanloc,'style','map');
    caxis([0 cup])
    
    subplot(2,4,fre+4)
    topoplot(mean(power_dz_rel,3),chanloc,'style','map');
    caxis([0 cup])
end


power_sum_tt = mean(spe_tt_rel,[1 3]);
power_sum_dz = mean(spe_dz_rel,[1 3]);
