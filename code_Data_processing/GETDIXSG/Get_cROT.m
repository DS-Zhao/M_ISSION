function [] = Get_cROT(path_obs,Sites_Info,sys)
%% Save the latitude, longitude, and cROT values of the puncture point
% INPUT:
%     path_obs: storage path of OBS files
%     Sites_Info: name and coordinate information of the stations
%     sys: Satellite system selection parameter status
%Save:
%     */ivcROT/syscROTdoy/sitedoysys_B_L_cROT.mat：latitude, longitude, cROT
%% written by Zhang P. et al., 2024/08    
%% --------------------------------------------------------------     
Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
current_doy=num2str(unique(doys));
list_obs=dir([path_obs,'/',current_doy,'/*.mat']);
path_ELE = fullfile(fileparts(path_obs), 'ivELE');
len=length(list_obs);
for i=1:len
    load([path_ELE,'/',current_doy,'/',list_obs(i).name],'-mat');%Load angle data
    load([path_obs,'/',current_doy,'/',list_obs(i).name],'-mat');%Load observation data
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
            [GPSB,GPSL,GPScROT]=GPS_cROT(obs,sx,sy,sz,EA.GPSel,EA.GPSaz);
            if all(all(GPScROT==0))
                continue;
            end
            GPS_B_L_cROT.B=GPSB;GPS_B_L_cROT.L=GPSL;GPS_B_L_cROT.cROT=GPScROT;
            gPath=fileparts(path_obs);
            if exist([gPath '/ivcROT' '/GPScROT' doy],'dir')==0
                mkdir([gPath '/ivcROT' '/GPScROT' doy]);
            end
            ffPath=[gPath '/ivcROT' '/GPScROT' doy '/' site doy 'GPS_B_L_cROT.mat'];
            save(ffPath,'GPS_B_L_cROT','-mat');
        end
    end
    
    if any(strcmp(fields, 'GLOC1P')) && any(strcmp(fields2, 'GLOel')) && sys(2)~=0
        if ~(all(all(obs.GLOC1P==0)) || all(all(obs.GLOC2P==0)))
            [GLOB,GLOL,GLOcROT]=GLO_cROT(obs,sx,sy,sz,EA.GLOel,EA.GLOaz);
            if all(all(GLOcROT==0))
                continue;
            end
            GLO_B_L_cROT.B=GLOB;GLO_B_L_cROT.L=GLOL;GLO_B_L_cROT.cROT=GLOcROT;
            gPath=fileparts(path_obs);
            if exist([gPath '/ivcROT' '/GLOcROT' doy],'dir')==0
                mkdir([gPath '/ivcROT' '/GLOcROT' doy]);
            end
            ffPath=[gPath '/ivcROT' '/GLOcROT' doy '/' site doy 'GLO_B_L_cROT.mat'];
            save(ffPath,'GLO_B_L_cROT','-mat');
        end
    end

     if any(strcmp(fields, 'BDSC2I')) && any(strcmp(fields2, 'BDSel')) && sys(4)~=0
        if ~(all(all(obs.BDSC2I==0)) || all(all(obs.BDSC6I==0)))
           [BDSB,BDSL,BDScROT]=BDS_cROT(obs,sx,sy,sz,EA.BDSel,EA.BDSaz);
            if all(all(BDScROT==0))
                continue;
            end
            BDS_B_L_cROT.B=BDSB;BDS_B_L_cROT.L=BDSL;BDS_B_L_cROT.cROT=BDScROT;
            gPath=fileparts(path_obs);
            if exist([gPath '/ivcROT' '/BDScROT' doy],'dir')==0
                mkdir([gPath '/ivcROT' '/BDScROT' doy]);
            end 
            ffPath=[gPath '/ivcROT' '/BDScROT' doy '/' site doy 'BDS_B_L_cROT.mat'];
            save(ffPath,'BDS_B_L_cROT','-mat');
        end
     elseif any(strcmp(fields, 'BDSC1I')) && any(strcmp(fields2, 'BDSel')) && sys(4)~=0
        if ~(all(all(obs.BDSC1I==0)) || all(all(obs.BDSC7I==0)))
            [BDSB,BDSL,BDScROT]=BDS1_cROT(obs,sx,sy,sz,EA.BDSel,EA.BDSaz);
            if all(all(BDScROT==0))
                continue;
            end
            BDS_B_L_cROT.B=BDSB;BDS_B_L_cROT.L=BDSL;BDS_B_L_cROT.cROT=BDScROT;
            gPath=fileparts(path_obs);
            if exist([gPath '/ivcROT' '/BDScROT' doy],'dir')==0   
                mkdir([gPath '/ivcROT' '/BDScROT' doy]);
            end
            ffPath=[gPath '/ivcROT' '/BDScROT' doy '/' site doy 'BDS_B_L_cROT.mat'];
            save(ffPath,'BDS_B_L_cROT','-mat');
        end
    end
    
    if any(strcmp(fields, 'GALC1C')) && any(strcmp(fields2, 'GALel')) && sys(3)~=0
        if ~(all(all(obs.GALC1C==0)) || all(all(obs.GALC5Q==0)))
            [GALB,GALL,GALcROT]=GAL_cROT(obs,sx,sy,sz,EA.GALel,EA.GALaz);
            if all(all(GALcROT==0))
                continue;
            end
            GAL_B_L_cROT.B=GALB;GAL_B_L_cROT.L=GALL;GAL_B_L_cROT.cROT=GALcROT;
            gPath=fileparts(path_obs);
            if exist([gPath '/ivcROT' '/GALcROT' doy],'dir')==0   
                mkdir([gPath '/ivcROT' '/GALcROT' doy]);
            end
            ffPath=[gPath '/ivcROT' '/GALcROT' doy '/' site doy 'GAL_B_L_cROT.mat'];
            save(ffPath,'GAL_B_L_cROT','-mat');
        end
    elseif any(strcmp(fields, 'GALC1X')) && any(strcmp(fields2, 'GALel')) && sys(3)~=0
        if ~(all(all(obs.GALC1X==0)) || all(all(obs.GALC5X==0)))
            [GALXB,GALXL,GALXcROT]=GALX_cROT(obs,sx,sy,sz,EA.GALel,EA.GALaz);
            if all(all(GALXcROT==0))
                continue;
            end
            GALcROT=GALXcROT;
            GALX_B_L_cROT.B=GALXB;GALX_B_L_cROT.L=GALXL;GAL_B_L_cROT.cROT=GALcROT;
            gPath=fileparts(path_obs);
            if exist([gPath '/ivcROT' '/GALcROT' doy],'dir')==0   
                mkdir([gPath '/ivcROT' '/GALcROT' doy]);
            end
            ffPath=[gPath '/ivcROT' '/GALcROT' doy '/' site doy 'GAL_B_L_cROT.mat'];
            save(ffPath,'GAL_B_L_cROT','-mat');
        end
    end
end
end
%% 
function [GPSB,GPSL,GPScROT]=GPS_cROT(obs,sx,sy,sz,E,A)
%% get GPSB,GPSL,GPScROT
% INPUT:
%      obs: struct of rinex files
%      sx,sy,sz:Receiver three coordinates
%      E:satellite elevation
%      A:satellite azimuth
% OUTPUT:
%     GPSB: latitude of Puncture point
%     GPSL: Longitude of puncture point
%% --------------------------------------------------------------------------
size2=size(obs.GPSC1W,2);
if size2<32
    obs.GPSL1C(:,size2+1:32)=0;obs.GPSL2W(:,size2+1:32)=0;
end
GPSL4=zeros(2880,size2);
c=299792458;                     %__________________________speed of light
GPS_f1=1575.42*10^6;                 %_________________________________unit:Hz
GPS_f2=1227.6*10^6;                  %_________________________________unit:Hz              
for i=1:size2 
    GPSL4(:,i)=(c/GPS_f1)*obs.GPSL1C(:,i)-(c/GPS_f2)*obs.GPSL2W(:,i);
    GPSL4(GPSL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(GPS_f1*GPS_f1*GPS_f2*GPS_f2)/((GPS_f1*GPS_f1-GPS_f2*GPS_f2)*Kei);%
    STEC=GPSL4*1e-3*para*1e+6;
    GPSdSTEC=nan(2880,size2);   
    for u=1:size(STEC,1)-1
        GPSdSTEC(u,:)=STEC(u+1,:)-STEC(u,:);
    end
end
[dif,GPSB,GPSL]=Get_dif(sx,sy,sz,E,A);

GPScROT=abs(GPSdSTEC./dif)*2;
end

%% ------------------------subfunction----------------------------
function [GLOB,GLOL,GLOcROT]=GLO_cROT(obs,sx,sy,sz,E,A)
%% get GLOB,GLOL,GLOcROT
% INPUT:
%      obs: struct of rinex files
%      sx,sy,sz:Receiver three coordinates
%      E:satellite elevation
%      A:satellite azimuth
% OUTPUT:
%     GLOB: latitude of Puncture point
%     GLOL: Longitude of puncture point
%% ------------------------------------------------------------------
size2=size(obs.GLOL1P,2);
if size2<24
    obs.GLOL1P(:,size2+1:24)=0;obs.GLOL2P(:,size2+1:24)=0;
end
GLOL4=zeros(2880,size2);
c=299792458;                     %__________________________speed of light
for i=1:size2
    Fre=[1,-4,5,6,1,-4, 5,6,-2,7,0,-1,-2,-7,0, -1, 4,-3,3, 2, 4, -3, 3, 2,0,0];
    GLO_f1(i)=(1602+Fre(i)*0.5625)*10^6;                 %_________________________________unit:Hz
    GLO_f2(i)=(1246+Fre(i)*0.4375)*10^6;                  %_________________________________unit:Hz
                   %_____________________wide lane ambiguity
end
for i=1:size2  %i is PRN number
    GLOL4(:,i)=(c/GLO_f1(i))*obs.GLOL1P(:,i)-(c/GLO_f2(i))*obs.GLOL2P(:,i);
    GLOL4(GLOL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(GLO_f1(i)*GLO_f1(i)*GLO_f2(i)*GLO_f2(i))/((GLO_f1(i)*GLO_f1(i)-GLO_f2(i)*GLO_f2(i))*Kei);
    GLOSTEC=GLOL4*para*1e+3;
    GLOdSTEC=nan(2880,size2);
    for u=1:size(GLOSTEC,1)-1
        GLOdSTEC(u,:)=GLOSTEC(u+1,:)-GLOSTEC(u,:);
    end
    [dif,GLOB,GLOL]=Get_dif(sx,sy,sz,E,A);
    GLOcROT=2*abs(GLOdSTEC./dif);
end
end

%% ---------------------------subfunction-------------------------
function [BDSB,BDSL,BDScROT]=BDS_cROT(obs,sx,sy,sz,E,A)
%% get BDSB,BDSL,BDScROT
% INPUT:
%      obs: struct of rinex files
%      sx,sy,sz:Receiver three coordinates
%      E:satellite elevation
%      A:satellite azimuth
% OUTPUT:
%     BDSB: latitude of Puncture point
%     BDSL: Longitude of puncture point
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
    BDSSTEC=BDSL4*para*1e+3;
    BDSdSTEC=nan(2880,size2);
    for u=1:size(BDSSTEC,1)-1
        BDSdSTEC(u,:)=BDSSTEC(u+1,:)-BDSSTEC(u,:);
    end
    [dif,BDSB,BDSL]=Get_dif(sx,sy,sz,E,A);
    BDScROT=2*abs(BDSdSTEC./dif);
end
end
%% ---------------------------subfunction-------------------------
function [BDSB,BDSL,BDScROT]=BDS1_cROT(obs,sx,sy,sz,E,A)
%% get BDSB,BDSL,BDScROT
% INPUT:
%      obs: struct of rinex files
%      sx,sy,sz:Receiver three coordinates
%      E:satellite elevation
%      A:satellite azimuth
% OUTPUT:
%     BDSB: latitude of Puncture point
%     BDSL: Longitude of puncture point
%% -----------------------------------------------------------------
    size2=size(obs.BDSL1I,2);
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
    BDSSTEC=BDSL4*para*1e+3;
    BDSdSTEC=nan(2880,size2);
    for u=1:size(BDSSTEC,1)-1
        BDSdSTEC(u,:)=BDSSTEC(u+1,:)-BDSSTEC(u,:);
    end
    [dif,BDSB,BDSL]=Get_dif(sx,sy,sz,E,A);
    BDScROT=2*abs(BDSdSTEC./dif);
end
end
%% --------------------------subfunction--------------------------
function [GALB,GALL,GALcROT]=GAL_cROT(obs,sx,sy,sz,E,A)
%% get GALB,GALL,GALcROT
% INPUT:
%      obs: struct of rinex files
%      sx,sy,sz:Receiver three coordinates
%      E:satellite elevation
%      A:satellite azimuth
% OUTPUT:
%     GALB: latitude of Puncture point
%     GALL: Longitude of puncture point
%% ------------------------------------------------------------------
size2=size(obs.GALC1C,2);
if size2<36
    obs.GALL1C(:,size2+1:36)=0;obs.GALL5Q(:,size2+1:36)=0;
end
GALL4=zeros(2880,size2);
c=299792458;                     %__________________________speed of light
GAL_f1=1575.42*10^6;                 %_________________________________unit:Hz
GAL_f5=1176.45*10^6;                  %_________________________________unit:Hz
for i=1:size2  %24  
    GALL4(:,i)=(c/GAL_f1)*obs.GALL1C(:,i)-(c/GAL_f5)*obs.GALL5Q(:,i);
    GALL4(GALL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(GAL_f1*GAL_f1*GAL_f5*GAL_f5)/((GAL_f1*GAL_f1-GAL_f5*GAL_f5)*Kei);
    GALSTEC=GALL4*para*1e+3;
    GALdSTEC=nan(2880,size2);
    for u=1:size(GALSTEC,1)-1
        GALdSTEC(u,:)=GALSTEC(u+1,:)-GALSTEC(u,:);
    end
    [dif,GALB,GALL]=Get_dif(sx,sy,sz,E,A);
    GALcROT=2*abs(GALdSTEC./dif);
end
end

%% --------------------------subfunction----------------------------
function [GALXB,GALXL,GALXcROT]=GALX_cROT(obs,sx,sy,sz,E,A)
%% get GALXB,GALXL,GALXcROT
% INPUT:
%      obs: struct of rinex files
%      sx,sy,sz:Receiver three coordinates
%      E:satellite elevation
%      A:satellite azimuth
% OUTPUT:
%     GALXB: latitude of Puncture point
%     GALXL: Longitude of puncture point
size2=size(obs.GALC1X,2);
if size2<36
    obs.GALL1X(:,size2+1:36)=0;obs.GALL5X(:,size2+1:36)=0;
end
GALXL4=zeros(2880,size2);
c=299792458;                     %__________________________speed of light
GALX_f1=1575.42*10^6;                 %_________________________________unit:Hz
GALX_f5=1176.45*10^6;                 %_________________________________unit:Hz
for i=1:size2  %i is PRN number
    GALXL4(:,i)=(c/GALX_f1)*obs.GALL1X(:,i)-(c/GALX_f5)*obs.GALL5X(:,i);
    GALXf1=1575.42*10^6;
    GALXf2=1176.45*10^6;
    GALXL4(GALXL4==0)=NaN;
    Kei=40.309*1e+16;
    para=(GALXf1*GALXf1*GALXf2*GALXf2)/((GALXf1*GALXf1-GALXf2*GALXf2)*Kei);
    GALXSTEC=GALXL4*para*1e+3;
    GALXdSTEC=nan(2880,size2);
    for u=1:size(GALXSTEC,1)-1
        GALXdSTEC(u,:)=GALXSTEC(u+1,:)-GALXSTEC(u,:);
    end
    [dif,GALXB,GALXL]=Get_dif(sx,sy,sz,E,A);
    GALXcROT=2*abs(GALXdSTEC./dif);
end
end

