clc;clear;

load('/Volumes/Backup Plus/疼痛数据/pname.mat');

fs = 512;

%% load the required geometrical information
        load('standard_bem.mat');
        load('elec_aligned.mat');
        
        headmodel_eeg = vol;
        
        headmodel_eeg = ft_convert_units(headmodel_eeg, 'mm');

filepath = '/Volumes/Backup Plus/疼痛数据/conn result new/';
pindex = 0;
Npain = 1;
Nhc = 1;

%% Start calculation
for p = 1:37
    if p < 20
        pname = ttname{p,1};
        if ismember(pname,{'CDX' 'TYL' 'ZQF' 'CZZ'})
            continue
        end
        pindex = pindex + 1;
        filename = [pname '_tt_conn_'];
        Npain = Npain + 1;
    else
        pname = dzname{p-19,1};
        pindex = pindex + 1;
        filename = [pname '_dz_conn_'];
        Nhc = Nhc + 1;
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
        
        % global efficiency
        glbeff_PLV(pindex,fre) = PLV.glbeff;
        glbeff_PLI(pindex,fre) = PLI.glbeff;
        glbeff_wPLI(pindex,fre) = wPLI.glbeff;
        
        % char path
        chpath_PLV(pindex,fre) = PLV.chpath;
        chpath_PLI(pindex,fre) = PLI.chpath;
        chpath_wPLI(pindex,fre) = wPLI.chpath;
        
        % small world
        smw_PLV(pindex,fre) = PLV.smwolrd;
        smw_PLI(pindex,fre) = PLI.smwolrd;
        smw_wPLI(pindex,fre) = wPLI.smwolrd;
        
        % global CC
        glbCC_PLV(pindex,fre) = PLV.glbCC;
        glbCC_PLI(pindex,fre) = PLI.glbCC;
        glbCC_wPLI(pindex,fre) = wPLI.glbCC;
        

    end
end


%% arrange mean and std
load('varanamelist.mat')



for n = 1: size(seqindex,1)
    paraname = seqindex{n,1};
    frename = seqindex{n,2};
    switch frename
        case 'delta'
            freindex = 1;
        case 'theta'
            freindex = 2;
        case 'alpha'
            freindex = 3;
        case 'beta'
            freindex = 4;
        case 'gamma'
            freindex = 5;
    end
    
    varname_PLV = [paraname '_PLV'];
    varname_PLI = [paraname '_PLI'];
    varname_wPLI = [paraname '_wPLI'];
    

    % pain
    var_temp = eval(varname_PLV);
    output_PLV_pain{n,1} = mean(var_temp(1:15,freindex),1);
    output_PLV_pain{n,2} = std(var_temp(1:15,freindex),1);
    output_PLV_pain{n,3} = [num2str(round(output_PLV_pain{n,1},3)) ' ± ' num2str(round(output_PLV_pain{n,2},3))];
    
    var_temp = eval(varname_PLI);
    output_PLI_pain{n,1} = mean(var_temp(1:15,freindex),1);
    output_PLI_pain{n,2} = std(var_temp(1:15,freindex),1);
    output_PLI_pain{n,3} = [num2str(round(output_PLI_pain{n,1},3)) ' ± ' num2str(round(output_PLI_pain{n,2},3))];
    
    var_temp = eval(varname_wPLI);
    output_wPLI_pain{n,1} = mean(var_temp(1:15,freindex),1);
    output_wPLI_pain{n,2} = std(var_temp(1:15,freindex),1);
    output_wPLI_pain{n,3} = [num2str(round(output_wPLI_pain{n,1},3)) ' ± ' num2str(round(output_wPLI_pain{n,2},3))];
    
    
    % control
    var_temp = eval(varname_PLV);
    output_PLV_hc{n,1} = mean(var_temp(16:33,freindex),1);
    output_PLV_hc{n,2} = std(var_temp(:,freindex),1);
    output_PLV_hc{n,3} = [num2str(round(output_PLV_hc{n,1},3)) ' ± ' num2str(round(output_PLV_hc{n,2},3))];
    
    var_temp = eval(varname_PLI);
    output_PLI_hc{n,1} = mean(var_temp(16:33,freindex),1);
    output_PLI_hc{n,2} = std(var_temp(:,freindex),1);
    output_PLI_hc{n,3} = [num2str(round(output_PLI_hc{n,1},3)) ' ± ' num2str(round(output_PLI_hc{n,2},3))];
    
    var_temp = eval(varname_wPLI);
    output_wPLI_hc{n,1} = mean(var_temp(16:33,freindex),1);
    output_wPLI_hc{n,2} = std(var_temp(:,freindex),1);
    output_wPLI_hc{n,3} = [num2str(round(output_wPLI_hc{n,1},3)) ' ± ' num2str(round(output_wPLI_hc{n,2},3))];
    
    
    
end
