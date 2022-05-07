
clc;clear;

load('/Volumes/Backup Plus/疼痛数据/pname.mat');

fs = 512;

%% load the required geometrical information
        load('standard_bem.mat');
        load('elec_aligned.mat');
        headmodel_eeg = vol;
        
        headmodel_eeg = ft_convert_units(headmodel_eeg, 'mm');

filepath = '/Volumes/Backup Plus/疼痛数据/conn result/';
savepath = '/Volumes/Backup Plus/疼痛数据/conn result new/';

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
        
        save([savepath filename frename '.mat'],'PLV','PLI','wPLI');
        
        clear data_complex  PLV PLI wPLI
        fprintf(['person ' num2str(p) 'fre ' num2str(fre) ' finished!\n'])
    end

end

        
               


%% Functions

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