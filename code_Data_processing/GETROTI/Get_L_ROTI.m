function [] = Get_L_ROTI(path_obs,Sites_Info,sys)
%% Input the OBS file and calculate ROTI based on carrier phase observation
% INPUT:
%     path_obs: storage path of OBS files
%     Sites_Info: name and coordinate information of the stations
%     sys: Satellite system selection parameter status
%Save:
%     */resROTI/sysROTIdoy/sitedoysysROTI.mat:ROIT calculation results
%% written by Zhao D,Zhang P. et al., 2024/08
%% ---------------------------------------------------------------------
Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
current_doy=num2str(unique(doys));
list_obs=dir([path_obs,'/',current_doy,'/*.mat']);
len=length(list_obs);

for i=1:len
    M_ROTI=[];
    load([path_obs,'/',current_doy,'/',list_obs(i).name],'-mat');
    site=list_obs(i).name(1:end-9);
    doy=list_obs(i).name(end-8:end-4);
    indices=doys==str2double(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);
    condition = all([sx == 0, sy == 0, sz == 0]);
    % If the condition is met, proceed to the next iteration
    if condition
        continue;
    end
    fields=fieldnames(obs);

    if any(strcmp(fields, 'GPSC1W' )) && sys(1)~=0
        if ~(all(all(obs.GPSC1W==0)) || all(all(obs.GPSC2W==0)))
            [GPSROTI]=GPS_ROTI(obs);
            figure;
            subplot(4,1,1);
            plot(obs.GPSL1C);
            subplot(4,1,2);
            plot(obs.GPSL2W);
            subplot(4,1,3);
            plot(GPSROTI);
            GPSROTI = checkAnomalies(GPSROTI);
            subplot(4,1,4);
            plot(GPSROTI);
            if all(all(GPSROTI==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resROTI' '/GPSROTI' doy],'dir')==0
                mkdir([gPath '/resROTI' '/GPSROTI' doy]);
            end
            ffPath=[gPath '/resROTI' '/GPSROTI' doy '/' site doy 'GPSROTI.mat'];
            save(ffPath,'GPSROTI','-mat');
            M_ROTI.GPSROTI=GPSROTI;
        end
    end

    if any(strcmp(fields, 'GLOC1P' ))&& sys(2)~=0
        if ~(all(all(obs.GLOC1P==0)) || all(all(obs.GLOC2P==0)))
            [GLOROTI]=GLO_ROTI(obs);
            GLOROTI = checkAnomalies(GLOROTI);
            if all(all(GLOROTI==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resROTI' '/GLOROTI' doy],'dir')==0
                mkdir([gPath '/resROTI' '/GLOROTI' doy]);
            end
            ffPath=[gPath '/resROTI' '/GLOROTI' doy '/' site doy 'GLOROTI.mat'];
            save(ffPath,'GLOROTI','-mat');
            M_ROTI.GLOROTI=GLOROTI;
        end
    end

    if any(strcmp(fields, 'BDSC2I' )) && sys(4)~=0
        if ~(all(all(obs.BDSC2I==0)) || all(all(obs.BDSC6I==0)))
            [BDSROTI]=BDS_ROTI(obs);
            BDSROTI = checkAnomalies(BDSROTI);
            if all(all(BDSROTI==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resROTI' '/BDSROTI' doy],'dir')==0
                mkdir([gPath '/resROTI' '/BDSROTI' doy]);
            end
            ffPath=[gPath '/resROTI' '/BDSROTI' doy '/' site doy 'BDSROTI.mat'];
            save(ffPath,'BDSROTI','-mat');
            M_ROTI.BDSROTI=BDSROTI;
        end
    elseif any(strcmp(fields, 'BDSC1I' )) && sys(4)~=0
        if ~(all(all(obs.BDSC1I==0)) || all(all(obs.BDSC7I==0)))
            [BDSROTI]=BDS1_ROTI(obs);
            BDSROTI = checkAnomalies(BDSROTI);
            if all(all(BDSROTI==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resROTI' '/BDSROTI' doy],'dir')==0
                mkdir([gPath '/resROTI' '/BDSROTI' doy]);
            end
            ffPath=[gPath '/resROTI' '/BDSROTI' doy '/' site doy 'BDSROTI.mat'];
            save(ffPath,'BDSROTI','-mat');
            M_ROTI.BDSROTI=BDSROTI;
        end
    end
    if any(strcmp(fields, 'GALC1C' )) && sys(3)~=0
        if ~(all(all(obs.GALC1C==0)) || all(all(obs.GALC5Q==0)))
            [GALROTI]=GAL_ROTI(obs);
            GALROTI = checkAnomalies(GALROTI);
            if all(all(GALROTI==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resROTI' '/GALROTI' doy],'dir')==0
                mkdir([gPath '/resROTI' '/GALROTI' doy]);
            end
            ffPath=[gPath '/resROTI' '/GALROTI' doy '/' site doy 'GALROTI.mat'];
            save(ffPath,'GALROTI','-mat');
            M_ROTI.GALROTI=GALROTI;
        end
    elseif any(strcmp(fields, 'GALC1X' ))  && sys(3)~=0
        if ~(all(all(obs.GALC1X==0)) || all(all(obs.GALC5X==0)))
            [GALROTI]=GALX_ROTI(obs);
            GALROTI = checkAnomalies(GALROTI);
            if all(all(GALROTI==0))
                continue;
            end
            gPath=fileparts(path_obs);
            if exist([gPath '/resROTI' '/GALROTI' doy],'dir')==0
                mkdir([gPath '/resROTI' '/GALROTI' doy]);
            end
            ffPath=[gPath '/resROTI' '/GALROTI' doy '/' site doy 'GALROTI.mat'];
            save(ffPath,'GALROTI','-mat');
            M_ROTI.GALXROTI=GALXROTI;
        end
    end
    if length(fieldnames(M_ROTI))>=2
    if exist([gPath '/resROTI' '/Multi_system_ROTI_in_one_mat'],'dir')==0
        mkdir([gPath '/resROTI' '/Multi_system_ROTI_in_one_mat']);
    end
    ffPath=[gPath '/resROTI' '/Multi_system_ROTI_in_one_mat' '/' site doy 'M_ROTI.mat'];
    save(ffPath,'M_ROTI','-mat');
    end

end
end
%% ---------------------------subfunction----------------------------------
function [ROTI]=GPS_ROTI(obs)
%% get  GPS ROTI observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     GPSROTI: ROTI calculation results
%% --------------------------------------------------------------------------
size2=size(obs.GPSL1C,2);
GPSL4=zeros(2880,size2);
ROTI=NaN(2880,size2);
c=299792458;                     %__________________________speed of light
GPS_f1=1575.42*10^6;                 %_________________________________unit:Hz
GPS_f2=1227.6*10^6;                  %_________________________________unit:Hz
for i=1:size2
    GPSL4(:,i)=(c/GPS_f1)*obs.GPSL1C(:,i)-(c/GPS_f2)*obs.GPSL2W(:,i);
    GPSL4(GPSL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(GPS_f1*GPS_f1*GPS_f2*GPS_f2)/((GPS_f1*GPS_f1-GPS_f2*GPS_f2)*Kei);
    STEC=GPSL4*para;
    dSTEC=nan(2880,size2);
    for u=1:size(STEC,1)-1
        dSTEC(u,:)=STEC(u+1,:)-STEC(u,:);
    end
    ROT=dSTEC*2;
    for m=11:2881
        ROTI(m-6,i)=std(ROT(m-10:m-1,i));
    end

end
end

%% ------------------------subfunction----------------------------
function [GLOROTI]=GLO_ROTI(obs)
%% get  GLO ROTI observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     GLOROTI: ROTI calculation results
%% --------------------------------------------------------------------------
size2=size(obs.GLOL1P,2);
if size2<24
    obs.GLOL1P(:,size2+1:24)=0;obs.GLOL2P(:,size2+1:24)=0;
end
GLOL4=zeros(2880,24);
GLOROTI=NaN(2880,24);
c=299792458;                     %__________________________speed of light
for i=1:size2
    Fre=[1,-4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2,0,0];
    GLO_f1(i)=(1602+Fre(i)*0.5625)*10^6;                 %_________________________________unit:Hz
    GLO_f2(i)=(1246+Fre(i)*0.4375)*10^6;                  %_________________________________unit:Hz
    %_____________________wide lane ambiguity
end
for i=1:24  %i is PRN number
    GLOL4(:,i)=(c/GLO_f1(i))*obs.GLOL1P(:,i)-(c/GLO_f2(i))*obs.GLOL2P(:,i);
    GLOL4(GLOL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(GLO_f1(i)*GLO_f1(i)*GLO_f2(i)*GLO_f2(i))/((GLO_f1(i)*GLO_f1(i)-GLO_f2(i)*GLO_f2(i))*Kei);
    GLOSTEC=GLOL4*para;
    GLOdSTEC=nan(2880,24);
    for u=1:size(GLOSTEC,1)-1
        GLOdSTEC(u,:)=GLOSTEC(u+1,:)-GLOSTEC(u,:);
    end
    GLOROT=GLOdSTEC*2;
    for m=11:2881
        GLOROTI(m-6,i)=std(GLOROT(m-10:m-1,i));
    end
end
end

%% ---------------------------subfunction-------------------------
function [BDSROTI]=BDS_ROTI(obs)
%% get  BDS ROTI observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     BDSROTI: ROTI calculation results
%% -----------------------------------------------------------------

size2=size(obs.BDSL2I,2);
BDS_f2=1561.098*10^6;
BDS_f6=1268.520*10^6;
BDSL4=zeros(2880,size2);
BDSROTI=NaN(2880,size2);
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
    BDSROT=BDSdSTEC*2;
    for m=11:2881
        BDSROTI(m-6,i)=std(BDSROT(m-10:m-1,i));
    end
end
end
%% ---------------------------subfunction-------------------------
function [BDSROTI]=BDS1_ROTI(obs)
%% get  BDS ROTI observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     BDSROTI: ROTI calculation results
%% -----------------------------------------------------------------
size2=size(obs.BDSC1I,2);
BDS_f1=1561.098*10^6;
BDS_f7=1207.140*10^6;
BDSL4=zeros(2880,size2);
BDSROTI=NaN(2880,size2);
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
    BDSROT=BDSdSTEC*2;
    for m=11:2881
        BDSROTI(m-6,i)=std(BDSROT(m-10:m-1,i));
    end
end
end
%% --------------------------subfunction--------------------------
function [GALROTI]=GAL_ROTI(obs)
%% get  GAL ROTI observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     GALROTI: ROTI calculation results
%% ------------------------------------------------------------------
size2=size(obs.GALL1C,2);
GALL4=zeros(2880,size2);
GALROTI=NaN(2880,size2);
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
    GALROT=GALdSTEC*2;
    for m=11:2881
        GALROTI(m-6,i)=std(GALROT(m-10:m-1,i));
    end
end
end

%% --------------------------subfunction----------------------------
function [GALXROTI]=GALX_ROTI(obs)
%% get  GALX ROTI observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     GALXROTI: ROTI calculation results
%% --------------------------------------------------------------------
size2=size(obs.GALL1X,2);
GALXL4=zeros(2880,size2);
GALXROTI=NaN(2880,size2);
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
    GALXROT=GALXdSTEC*2;
    for m=11:2881
        GALXROTI(m-6,i)=std(GALXROT(m-10:m-1,i));
    end
end
end


