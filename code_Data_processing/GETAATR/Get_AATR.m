function [] = Get_AATR(path_obs,Sites_Info,sys)
%% Input the OBS file and calculate AATR based on carrier phase observation
% INPUT:
%     path_obs: storage path of OBS files
%     Sites_Info: name and coordinate information of the stations
%     sys: Satellite system selection parameter status
%Save:
%     */resAATR/sysAATRdoy/sitedoysysAATR.mat:AATR calculation results
%     */resAATR/Multi_system_AATR_in_one_mat/sitedoyM_AATR.mat: Multi system RMS aAATR calculation results
%     */resRMSAATR/sysRMSAATRdoy/sitedoysysRMSAATR.mat: RMS aAATR calculation results
%     */resRMSAATR/Multi_system_RMSAATR_in_one_mat/sitedoyM_RAATR.mat: Multi system RMS aAATR calculation results
%% written by Zhang X,Zhang P. et al., 2024/08
%% ---------------------------------------------------------------------
Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
current_doy=num2str(unique(doys));
list_obs=dir([path_obs,'/',current_doy,'/*.mat']);
path_ELE = fullfile(fileparts(path_obs), 'ivELE');
len=length(list_obs);

for i=1:len
    M_RAATR=[];
    load([path_ELE,'/',current_doy,'/',list_obs(i).name],'-mat');
    load([path_obs,'/',current_doy,'/',list_obs(i).name],'-mat');
    site=list_obs(i).name(1:end-9);
    doy=list_obs(i).name(end-8:end-4);
    indices=doys==str2double(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);
    condition = all([sx == 0, sy == 0, sz == 0]);
    if condition
        continue;
    end
    fields=fieldnames(obs);
    fields2=fieldnames(EA);

    if any(strcmp(fields, 'GPSC1W')) && any(strcmp(fields2, 'GPSel')) && sys(1)~=0
        if ~(all(all(obs.GPSC1W==0)) || all(all(obs.GPSC2W==0)))
            [GPSAATR,GPSRAATR]=GPS_AATR(obs,EA.GPSel);
            if all(all(GPSRAATR==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resAATR' '/GPSAATR' doy],'dir')==0
                mkdir([gPath '/resAATR' '/GPSAATR' doy]);
            end
            ffPath=[gPath '/resAATR' '/GPSAATR' doy '/' site doy 'GPSAATR.mat'];
            save(ffPath,'GPSAATR','-mat');
            M_AATR.GPSAATR=GPSAATR;
            if exist([gPath '/resRMSAATR' '/GPSRMSAATR' doy],'dir')==0
                mkdir([gPath '/resRMSAATR' '/GPSRMSAATR' doy]);
            end
            ffPath=[gPath '/resRMSAATR' '/GPSRMSAATR' doy '/' site doy 'GPSRMSAATR.mat'];
            save(ffPath,'GPSRAATR','-mat');
            M_RAATR.GPSRAATR=GPSRAATR;
        end
    end

    if any(strcmp(fields, 'GLOC1P')) && any(strcmp(fields2, 'GLOel')) && sys(2)~=0
        if ~(all(all(obs.GLOC1P==0)) || all(all(obs.GLOC2P==0)))
            [GLOAATR,GLORAATR]=GLO_AATR(obs,EA.GLOel);
            if all(all(GLORAATR==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resAATR' '/GLOAATR' doy],'dir')==0
                mkdir([gPath '/resAATR' '/GLOAATR' doy]);
            end
            ffPath=[gPath '/resAATR' '/GLOAATR' doy '/' site doy 'GLOAATR.mat'];
            save(ffPath,'GLOAATR','-mat');
            M_AATR.GLOAATR=GLOAATR;
            if exist([gPath '/resRMSAATR' '/GLORMSAATR' doy],'dir')==0
                mkdir([gPath '/resRMSAATR' '/GLORMSAATR' doy]);
            end
            ffPath=[gPath '/resRMSAATR' '/GLORMSAATR' doy '/' site doy 'GLORMSAATR.mat'];
            save(ffPath,'GLORAATR','-mat');
            M_RAATR.GLORAATR=GLORAATR;
        end
    end

    if any(strcmp(fields, 'BDSC2I')) && any(strcmp(fields2, 'BDSel')) && sys(4)~=0
        if ~(all(all(obs.BDSC2I==0)) || all(all(obs.BDSC6I==0)))
            [BDSAATR,BDSRAATR]=BDS_AATR(obs,EA.BDSel);
            if all(all(BDSRAATR==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resAATR' '/BDSAATR' doy],'dir')==0
                mkdir([gPath '/resAATR' '/BDSAATR' doy]);
            end
            ffPath=[gPath '/resAATR' '/BDSAATR' doy '/' site doy 'BDSAATR.mat'];
            save(ffPath,'BDSAATR','-mat');
            M_AATR.BDSAATR=BDSAATR;
            if exist([gPath '/resRMSAATR' '/BDSRMSAATR' doy],'dir')==0
                mkdir([gPath '/resRMSAATR' '/BDSRMSAATR' doy]);
            end
            ffPath=[gPath '/resRMSAATR' '/BDSRMSAATR' doy '/' site doy 'BDSRMSAATR.mat'];
            save(ffPath,'BDSRAATR','-mat');
            M_RAATR.BDSRAATR=BDSRAATR;
        end
    elseif any(strcmp(fields, 'BDSC1I')) && any(strcmp(fields2, 'BDSel')) && sys(4)~=0
        if ~(all(all(obs.BDSC1I==0)) || all(all(obs.BDSC7I==0)))
            [BDSAATR,BDSRAATR]=BDS1_AATR(obs,EA.BDSel);
            if all(all(BDSRAATR==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resAATR' '/BDSAATR' doy],'dir')==0
                mkdir([gPath '/resAATR' '/BDSAATR' doy]);
            end
            ffPath=[gPath '/resAATR' '/BDSAATR' doy '/' site doy 'BDSAATR.mat'];
            save(ffPath,'BDSAATR','-mat');
            M_AATR.BDSAATR=BDSAATR;
            if exist([gPath '/resRMSAATR' '/BDSRMSAATR' doy],'dir')==0
                mkdir([gPath '/resRMSAATR' '/BDSRMSAATR' doy]);
            end
            ffPath=[gPath '/resRMSAATR' '/BDSRMSAATR' doy '/' site doy 'BDSRMSAATR.mat'];
            save(ffPath,'BDSRAATR','-mat');
            M_RAATR.BDSRAATR=BDSRAATR;
        end
    end

    if any(strcmp(fields, 'GALC1C')) && any(strcmp(fields2, 'GALel')) && sys(3)~=0
        if ~(all(all(obs.GALC1C==0)) || all(all(obs.GALC5Q==0)))
            [GALAATR,GALRAATR]=GAL_AATR(obs,EA.GALel);
            if all(all(GALRAATR==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resAATR' '/GALAATR' doy],'dir')==0
                mkdir([gPath '/resAATR' '/GALAATR' doy]);
            end
            ffPath=[gPath '/resAATR' '/GALAATR' doy '/' site doy 'GALAATR.mat'];
            save(ffPath,'GALAATR','-mat');
            M_AATR.GALAATR=GALAATR;
            if exist([gPath '/resRMSAATR' '/GALRMSAATR' doy],'dir')==0
                mkdir([gPath '/resRMSAATR' '/GALRMSAATR' doy]);
            end
            ffPath=[gPath '/resRMSAATR' '/GALRMSAATR' doy '/' site doy 'GALRMSAATR.mat'];
            save(ffPath,'GALRAATR','-mat');
            M_RAATR.GALRAATR=GALRAATR;
        end
    elseif any(strcmp(fields, 'GALC1X')) && any(strcmp(fields2, 'GALel')) && sys(3)~=0
        if ~(all(all(obs.GALC1X==0)) || all(all(obs.GALC5X==0)))
            [GALXAATR,GALXRAATR]=GALX_AATR(obs,EA.GALel);
            if all(all(GALXRAATR==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resAATR' '/GALXAATR' doy],'dir')==0
                mkdir([gPath '/resAATR' '/GALXAATR' doy]);
            end
            ffPath=[gPath '/resAATR' '/GALXAATR' doy '/' site doy 'GALXAATR.mat'];
            save(ffPath,'GALXAATR','-mat');
            M_AATR.GALXAATR=GALXAATR;
            if exist([gPath '/resRMSAATR' '/GALXRMSAATR' doy],'dir')==0
                mkdir([gPath '/resRMSAATR' '/GALXRMSAATR' doy]);
            end
            ffPath=[gPath '/resRMSAATR' '/GALXRMSAATR' doy '/' site doy 'GALXRMSAATR.mat'];
            save(ffPath,'GALXRAATR','-mat');
            M_RAATR.GALXRAATR=GALXRAATR;
        end
    end
    if length(fieldnames(M_AATR))>=2
        if exist([gPath '/resAATR' '/Multi_system_AATR_in_one_mat'],'dir')==0
            mkdir([gPath '/resAATR' '/Multi_system_AATR_in_one_mat']);
        end
        ffPath=[gPath '/resAATR' '/Multi_system_AATR_in_one_mat' '/' site doy 'M_AATR.mat'];
        save(ffPath,'M_AATR','-mat');
    end
    if length(fieldnames(M_RAATR))>=2
        if exist([gPath '/resRMSAATR' '/Multi_system_RMSAATR_in_one_mat'],'dir')==0
            mkdir([gPath '/resRMSAATR' '/Multi_system_RMSAATR_in_one_mat']);
        end
        ffPath=[gPath '/resRMSAATR' '/Multi_system_RMSAATR_in_one_mat' '/' site doy 'M_RAATR.mat'];
        save(ffPath,'M_RAATR','-mat');
    end
end
end
%% ---------------------------subfunction----------------------------------
function [GPSAATR,GPSRAATR]=GPS_AATR(obs,EL)
%% get  GPS AATR observations
% INPUT:
%      obs: struct of rinex files
%      EL:satellite azimuth and elevation
% OUTPUT:
%     GPSAATR: AATR calculation results
%     GPSRAATR:RMS AATR calculation results
%% --------------------------------------------------------------------------
size2=size(obs.GPSL1C,2);
GPSL4=zeros(2880,size2);
c=299792458;                     %__________________________speed of light
GPS_f1=1575.42*10^6;                 %_________________________________unit:Hz
GPS_f2=1227.6*10^6;                  %_________________________________unit:Hz
for i=1:size2
    GPSL4(:,i)=(c/GPS_f1)*obs.GPSL1C(:,i)-(c/GPS_f2)*obs.GPSL2W(:,i);
    GPSL4(GPSL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(GPS_f1*GPS_f1*GPS_f2*GPS_f2)/((GPS_f1*GPS_f1-GPS_f2*GPS_f2)*Kei);
    STEC=GPSL4*para;
    GPSdSTEC=nan(2880,size2);
    for u=1:size(STEC,1)-1
        GPSdSTEC(u,:)=STEC(u+1,:)-STEC(u,:);
    end
    el=EL;
    el(el==0)=nan;
    M=1./(1-(6371/(350+6371))^2*(cos(el).*cos(el)));
    GPSAATR=2*GPSdSTEC./M;
    GPSAATR = checkAnomalies(GPSAATR);
    GPSRAATR=zeros(24,1);
    j=1;
    for i=121:120:2881
        GPSRAATR(j)=rms(GPSAATR(i-120:i-1,:),"all",'omitnan');
        j=j+1;
    end
end
end

%% ------------------------subfunction----------------------------
function [GLOAATR,GLORAATR]=GLO_AATR(obs,EL)
%% get  GLO AATR observations
% INPUT:
%      obs: struct of rinex files
%      EL:satellite azimuth and elevation
% OUTPUT:
%     GLOAATR: AATR calculation results
%     GLORAATR:RMS AATR calculation results
%% ------------------------------------------------------------------
size2=size(obs.GLOL1P,2);
GLOL4=zeros(2880,size2);
c=299792458;                     %__________________________speed of light
for i=1:size2
    Fre=[1,-4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2,0,0];
    %      Fre=[1,-4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2];
    GLO_f1(i)=(1602+Fre(i)*0.5625)*10^6;                 %_________________________________unit:Hz
    GLO_f2(i)=(1246+Fre(i)*0.4375)*10^6;                  %_________________________________unit:Hz
    %_____________________wide lane ambiguity
end
for i=1:size2  %i is PRN number
    GLOL4(:,i)=(c/GLO_f1(i))*obs.GLOL1P(:,i)-(c/GLO_f2(i))*obs.GLOL2P(:,i);
    GLOL4(GLOL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(GLO_f1(i)*GLO_f1(i)*GLO_f2(i)*GLO_f2(i))/((GLO_f1(i)*GLO_f1(i)-GLO_f2(i)*GLO_f2(i))*Kei);
    GLOSTEC=GLOL4*para;
    GLOdSTEC=nan(2880,size2);
    for u=1:size(GLOSTEC,1)-1
        GLOdSTEC(u,:)=GLOSTEC(u+1,:)-GLOSTEC(u,:);
    end
    el=EL;
    el(el==0)=nan;
    M=1./(1-(6371/(350+6371))^2*(cos(el).*cos(el)));
    GLOAATR=2*GLOdSTEC./M;
    GLORAATR=zeros(24,1);
    j=1;
    for i=121:120:2881
        GLORAATR(j)=rms(GLOAATR(i-120:i-1,:),"all",'omitnan');
        j=j+1;
    end
end
end

%% ---------------------------subfunction-------------------------
function [BDSAATR,BDSRAATR]=BDS_AATR(obs,EL)
%% get  BDS AATR observations
% INPUT:
%      obs: struct of rinex files
%      EL:satellite azimuth and elevation
% OUTPUT:
%     BDSAATR: AATR calculation results
%     BDSRAATR:RMS AATR calculation results
%% -----------------------------------------------------------------
size2=size(obs.BDSL2I,2);
BDS_f2=1561.098*10^6;
BDS_f6=1268.520*10^6;
BDSL4=zeros(2880,size2);
c=299792458;                     %__________________________speed of light
for i=1:size2  %i is PRN number
    BDSL4(:,i)=(c/BDS_f2)*obs.BDSL2I(:,i)-(c/BDS_f6)*obs.BDSL6I(:,i);
    BDSL4(BDSL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(BDS_f2*BDS_f2*BDS_f6*BDS_f6)/((BDS_f2*BDS_f2-BDS_f6*BDS_f6)*Kei);
    BDSSTEC=BDSL4*para;
    BDSdSTEC=nan(2880,size2);
    for u=1:size(BDSSTEC,1)-1
        BDSdSTEC(u,:)=BDSSTEC(u+1,:)-BDSSTEC(u,:);
    end
    el=EL;
    el(el==0)=nan;
    M=1./(1-(6371/(350+6371))^2*(cos(el).*cos(el)));
    BDSAATR=2*BDSdSTEC./M;
    BDSRAATR=zeros(24,1);
    j=1;
    for i=121:120:2881
        BDSRAATR(j)=rms(BDSAATR(i-120:i-1,:),"all",'omitnan');
        j=j+1;
    end
end
end
%% ---------------------------subfunction-------------------------
function [BDSAATR,BDSRAATR]=BDS1_AATR(obs,EL)
%% get  BDS AATR observations
% INPUT:
%      obs: struct of rinex files
%      EL:satellite azimuth and elevation
% OUTPUT:
%     BDSAATR: AATR calculation results
%     BDSRAATR:RMS AATR calculation results
%% -----------------------------------------------------------------
size2=size(obs.BDSC1I,2);
BDS_f1=1561.098*10^6;
BDS_f7=1207.140*10^6;
if size2<61
    obs.BDSL1I(:,size2+1:61)=0;
    obs.BDSL7I(:,size2+1:61)=0;
end
BDSL4=zeros(2880,size2);
c=299792458;                     %__________________________speed of light
for i=1:size2  %i is PRN number
    BDSL4(:,i)=(c/BDS_f1)*obs.BDSL1I(:,i)-(c/BDS_f7)*obs.BDSL7I(:,i);
    BDSL4(BDSL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(BDS_f1*BDS_f1*BDS_f7*BDS_f7)/((BDS_f1*BDS_f1-BDS_f7*BDS_f7)*Kei);
    BDSSTEC=BDSL4*para;
    BDSdSTEC=nan(2880,size2);
    for u=1:size(BDSSTEC,1)-1
        BDSdSTEC(u,:)=BDSSTEC(u+1,:)-BDSSTEC(u,:);
    end
    el=EL;
    el(el==0)=nan;
    M=1./(1-(6371/(350+6371))^2*(cos(el).*cos(el)));
    BDSAATR=2*BDSdSTEC./M;
    BDSRAATR=zeros(24,1);
    j=1;
    for i=121:120:2881
        BDSRAATR(j)=rms(BDSAATR(i-120:i-1,:),"all",'omitnan');
        j=j+1;
    end
end
end
%% --------------------------subfunction--------------------------
function [GALAATR,GALRAATR]=GAL_AATR(obs,EL)
%% get  GAL AATR observations
% INPUT:
%      obs: struct of rinex files
%      EL:satellite azimuth and elevation
% OUTPUT:
%     GALAATR: AATR calculation results
%     GALRAATR:RMS AATR calculation results
%% ------------------------------------------------------------------
size2=size(obs.GALL1C,2);
GALL4=zeros(2880,size2);
c=299792458;                     %__________________________speed of light
GAL_f1=1575.42*10^6;                 %_________________________________unit:Hz
GAL_f5=1176.45*10^6;                  %_________________________________unit:Hz
for i=1:size2  %24
    GALL4(:,i)=(c/GAL_f1)*obs.GALL1C(:,i)-(c/GAL_f5)*obs.GALL5Q(:,i);
    GALL4(GALL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(GAL_f1*GAL_f1*GAL_f5*GAL_f5)/((GAL_f1*GAL_f1-GAL_f5*GAL_f5)*Kei);
    GALSTEC=GALL4*para;
    GALdSTEC=nan(2880,size2);
    for u=1:size(GALSTEC,1)-1
        GALdSTEC(u,:)=GALSTEC(u+1,:)-GALSTEC(u,:);
    end
    el=EL;
    el(el==0)=nan;
    M=1./(1-(6371/(350+6371))^2*(cos(el).*cos(el)));
    GALAATR=2*GALdSTEC./M;
    GALRAATR=zeros(24,1);
    j=1;
    for i=121:120:2881
        GALRAATR(j)=rms(GALAATR(i-120:i-1,:),"all",'omitnan');
        j=j+1;
    end
end
end

%% --------------------------subfunction----------------------------
function [GALXAATR,GALXRAATR]=GALX_AATR(obs,EL)
%% get  GALX AATR observations
% INPUT:
%      obs: struct of rinex files
%      EL:satellite azimuth and elevation
% OUTPUT:
%     GALXAATR: AATR calculation results
%     GALXRAATR:RMS AATR calculation results
%% --------------------------------------------------------------------
size2=size(obs.GALL1X,2);
GALXL4=zeros(2880,size2);
c=299792458;                     %__________________________speed of light
GALX_f1=1575.42*10^6;                 %_________________________________unit:Hz
GALX_f2=1176.45*10^6;                 %_________________________________unit:Hz
for i=1:size2  %i is PRN number
    GALXL4(:,i)=(c/GALX_f1)*obs.GALL1X(:,i)-(c/GALX_f2)*obs.GALL5X(:,i);
    GALXL4(GALXL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(GALX_f1*GALX_f1*GALX_f2*GALX_f2)/((GALX_f1*GALX_f1-GALX_f2*GALX_f2)*Kei);
    GALXSTEC=GALXL4*para;
    GALXdSTEC=nan(2880,size2);
    for u=1:size(GALXSTEC,1)-1
        GALXdSTEC(u,:)=GALXSTEC(u+1,:)-GALXSTEC(u,:);
    end
    el=EL;
    el(el==0)=nan;
    M=1./(1-(6371/(350+6371))^2*(cos(el).*cos(el)));
    GALXAATR=2*GALXdSTEC./M;
    GALXRAATR=zeros(24,1);
    j=1;
    for i=121:120:2881
        GALXRAATR(j)=rms(GALXAATR(i-120:i-1,:),"all",'omitnan');
        j=j+1;
    end
end
end
