clc;clear;

load('/Volumes/Backup Plus/疼痛数据/pname.mat');

fs = 512;

%% load the required geometrical information

filepath = '/Volumes/Backup Plus/疼痛数据/conn result new/';

%% Start calculation
for p = 1:37
    if p < 20
        pname = ttname{p,1};
        if ismember(pname,{'CDX' 'TYL' 'ZQF' 'CZZ'})
            continue
        end
        filename = [pname '_tt_conn_'];
    else
        pname = dzname{p-19,1};
        filename = [pname '_dz_conn_'];
    end
    
    
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
        load([filepath filename frename '.mat']);  % parameter: source
        
        %%

                
        % charpath
        plvD_temp = distance_bin(PLV.conn_thre_b);
        PLV.chpath = charpath(plvD_temp,0,0);
        
        pliD_temp = distance_bin(PLI.conn_thre_b);
        PLI.chpath = charpath(pliD_temp,0,0);
        
        wpliD_temp = distance_bin(wPLI.conn_thre_b);
        wPLI.chpath = charpath(wpliD_temp,0,0);
        

        
        clear data_complex  PLV PLI wPLI
        fprintf(['person ' num2str(p) 'fre ' num2str(fre) ' finished!\n'])
    end

end