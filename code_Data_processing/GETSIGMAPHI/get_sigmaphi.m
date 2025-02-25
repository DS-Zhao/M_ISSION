function get_sigmaphi(filtered_data,options,folderPath)
%% Calculate sigmaphi
% INPUT:
%     filtered_data: Filtered data
%     options: PPPH setting parameters
%     folderPath：folder path
% SAVE:
%     */resSIGMAPHI/syssigmaphidoy/sitedoysyssigmaphi.mat:sigmaphi calculation results
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
if options.system.gps==1 && (options.wave(1)~=0||options.wave(2)~=0)
    if options.wave(1)==1 && any(any(~isnan(filtered_data.obs.l1(:,1:32))))
        GPSdata=filtered_data.obs.l1(:,1:32);
        GPSsigmaphi.L1=sigmaphi(GPSdata);
    end
    if options.wave(2)==1 && any(any(~isnan(filtered_data.obs.l2(:,1:32))))
        GPSdata=filtered_data.obs.l2(:,1:32);
        GPSsigmaphi.L2=sigmaphi(GPSdata);
    end
    if exist("GPSsigmaphi")
    grandPath=fileparts(fileparts(folderPath));
    year=num2str(filtered_data.inf.time.first(1));
    doy=filtered_data.inf.time.doy;
    [~, fileName, ~] = fileparts(filtered_data.files.rinex);
    name=[fileName(1:4),year(3:4),num2str(doy,'%03d'),'GPSsigmaphi.mat'];
    newDirectoryName = ['resSIGMAPHI/','GPSsigmaphi',year(3:4),num2str(doy,'%03d')];
    cordata_path= fullfile(grandPath, newDirectoryName);
    if exist(cordata_path,'dir')==0
        mkdir(cordata_path);
    end
    save([cordata_path,'\',name],'GPSsigmaphi','-mat');
end
end
if options.system.glo==1 && (options.wave(1)~=0||options.wave(2)~=0)
    if options.wave(1)==1 && any(any(~isnan(filtered_data.obs.l1(:,33:58))))
        GLOdata=filtered_data.obs.l1(:,33:58);
        GLOsigmaphi.L1=sigmaphi(GLOdata);
    end
    if options.wave(2)==1 && any(any(~isnan(filtered_data.obs.l2(:,33:58))))
        GLOdata=filtered_data.obs.l2(:,33:58);
        GLOsigmaphi.L2=sigmaphi(GLOdata);
    end
    if exist("GLOsigmaphi")
    grandPath=fileparts(fileparts(folderPath));
    year=num2str(filtered_data.inf.time.first(1));
    doy=filtered_data.inf.time.doy;
    [~, fileName, ~] = fileparts(filtered_data.files.rinex);
    name=[fileName(1:4),year(3:4),num2str(doy,'%03d'),'GLOsigmaphi.mat'];
    newDirectoryName = ['resSIGMAPHI/','GLOsigmaphi',year(3:4),num2str(doy,'%03d')];
    cordata_path= fullfile(grandPath, newDirectoryName);
    if exist(cordata_path,'dir')==0
        mkdir(cordata_path);
    end
    save([cordata_path,'\',name],'GLOsigmaphi','-mat');
    end
end
end
function sigmaphi=sigmaphi(data)
[~,q]=size(data);
sigmaphi=nan(2880,q);
for i=1:q
    for m=11:2881
        sigmaphi(m-6,i)=std(data(m-10:m-1,i),1);
    end
end
end
