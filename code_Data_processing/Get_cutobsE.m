function [] = Get_cutobsE(path_obs,path_sp3,Sites_Info,lim,sate_mark)
%%Delete data below the cutoff angle and standardize the size of the observation data array 
%% and output angle data
% INPUT:
%     path_obs: storage path of obs files
%     path_sp3: storage path of sp3 files
%     Sites_Info: name and coordinate information of the stations
%     lim: cut-off angle
%     sate_mark: satellite status identification
%Save:
%     */OBS_cut/doy/sitedoy.mat：Processed OBS data
%     */ivELE/doy/sitedoy.mat:EA.el&&EA.az
%% written by Jin R et al., 2013/5/10, doi:10.1007/s10291-012-0279-3
%% modified by Zhou C. et al., 2021/12/12
%% modified by Zhang P. et al., 2024/08
%% ---------------------------------------------------------------------
Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
current_doy=num2str(unique(doys));
list_obs=dir([path_obs,'/',current_doy,'/*.mat']);
list_sp3=dir([path_sp3 '/*.mat']);
len=length(list_obs);
for i=1:len
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

    sp3index=find_sp3(list_sp3,doy);
    load([path_sp3 '/' list_sp3(sp3index).name],'-mat');
    fields=fieldnames(obs);
    fields2=fieldnames(sate);
    if any(~isnan(find(strcmp(fields, 'GPSL1C')))) && any(~isnan(find(strcmp(fields2, 'gpsx'))))
        gps=size(obs.GPSL1C,2);
        if gps<=size(sate.gpsx,2)
        sate.gpsx(:,gps+1:end)=[]; sate.gpsy(:,gps+1:end)=[]; sate.gpsz(:,gps+1:end)=[];
        end
    end

    if any(~isnan(find(strcmp(fields, 'GLOL1P' )))) && any(~isnan(find(strcmp(fields2, 'glox' ))))
        glo=size(obs.GLOL1P,2);
        if glo<=size(sate.glox,2)
        sate.glox(:,glo+1:end)=[]; sate.gloy(:,glo+1:end)=[]; sate.gloz(:,glo+1:end)=[];
        end
    end
    if any(~isnan(find(strcmp(fields, 'GALL1C' )))) && any(~isnan(find(strcmp(fields2, 'galx' ))))
        gal=size(obs.GALL1C,2);
        if gal<=size(sate.galx,2)
        sate.galx(:,gal+1:end)=[]; sate.galy(:,gal+1:end)=[]; sate.galz(:,gal+1:end)=[];
        end
    elseif any(~isnan(find(strcmp(fields, 'GALL1X' )))) && any(~isnan(find(strcmp(fields2, 'galx' ))))
        gal=size(obs.GALL1X,2);
        if gal<=size(sate.galx,2)
        sate.galx(:,gal+1:end)=[]; sate.galy(:,gal+1:end)=[]; sate.galz(:,gal+1:end)=[];
        end
    end

    if any(~isnan(find(strcmp(fields, 'BDSL2I' )))) && any(~isnan(find(strcmp(fields2, 'bdsx' ))))
        bds=size(obs.BDSL2I,2);
        if bds<=size(sate.bdsx,2)
        sate.bdsx(:,bds+1:end)=[]; sate.bdsy(:,bds+1:end)=[]; sate.bdsz(:,bds+1:end)=[];
        end
    elseif any(~isnan(find(strcmp(fields, 'BDSL1I' )))) && any(~isnan(find(strcmp(fields2, 'bdsx' ))))
        bds=size(obs.BDSL1I,2);
        if bds<=size(sate.bdsx,2)
        sate.bdsx(:,bds+1:end)=[]; sate.bdsy(:,bds+1:end)=[]; sate.bdsz(:,bds+1:end)=[];
        end
    end
    [EA,obs]=cut_obs(sate,sx,sy,sz,obs,lim,sate_mark);
    path_cutobs = fullfile(fileparts(path_obs), 'raw_OBS_cut',doy); 
    if exist(path_cutobs,'dir')==0
    mkdir(path_cutobs); 
    end
    filename2 = fullfile(path_cutobs, list_obs(i).name);
    coor=Coor(index,1:3);
    save(filename2, 'obs','coor','EA', '-mat');
end
end
%% ----------------subfunction-----------------
function [EA,obs]=cut_obs(sate,sx,sy,sz,obs,lim,sate_mark)
%% get the observations under specified cut-off angle and output angle data
% INPUT: 
%     sate: precise coordinates of the satellites
%     sx: X coordinate of the station
%     sy: Y coordinate of the station
%     sz: Z coordinate of the station
%     obs: original observation structs
%     lim: cut-off angle
%     sate_mark: satellite status identification
% OUTPUT：
%      obs: updated observation structs
%      EA:satellite azimuth and elevation
fields=fieldnames(obs);
% cut gps obs
if ~isnan(find(strcmp(fields, 'GPSC1W' )))
    gpsline=size(obs.GPSC2W,1);
    gpssate=size(sate_mark.gps,2);
    if gpsline<2880
        obs.GPSL1C(gpsline+1:2880,:)=0;obs.GPSL2W(gpsline+1:2880,:)=0;
        obs.GPSC1W(gpsline+1:2880,:)=0;obs.GPSC2W(gpsline+1:2880,:)=0;
    end
    gpsl=size(obs.GPSC2W,2);
    if gpsl<gpssate
        obs.GPSL1C(:,gpsl+1:gpssate)=0;obs.GPSL2W(:,gpsl+1:gpssate)=0;
        obs.GPSC1W(:,gpsl+1:gpssate)=0;obs.GPSC2W(:,gpsl+1:gpssate)=0;
    end
    
    gpsdelete=find(sate_mark.gps==0);
    gpsnum=size(sate.gpsx,2);    
    if gpsnum<size(obs.GPSC1W,2)
        if ~isnan(gpsdelete)
            k=length(gpsdelete);
            for i=k:-1:1
                obs.GPSC1W(:,gpsdelete(i))=0;
                obs.GPSC2W(:,gpsdelete(i))=0;
                obs.GPSL1C(:,gpsdelete(i))=0;
                obs.GPSL2W(:,gpsdelete(i))=0;
            end
        end
        obs.GPSC1W=obs.GPSC1W(:,1:gpsnum);
        obs.GPSC2W=obs.GPSC2W(:,1:gpsnum);
        obs.GPSL1C=obs.GPSL1C(:,1:gpsnum);
        obs.GPSL2W=obs.GPSL2W(:,1:gpsnum);
    end
    obs.GPSC1W(isnan(obs.GPSC1W))=0;
    obs.GPSC2W(isnan(obs.GPSC2W))=0;
    obs.GPSL1C(isnan(obs.GPSL1C))=0;
    obs.GPSL2W(isnan(obs.GPSL2W))=0;
    
    gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;

    for i=1:gpsnum
        for j=1:2880
            if obs.GPSL1C(j,i)==0 || obs.GPSL2W(j,i)==0 || obs.GPSC1W(j,i)==0 || obs.GPSC2W(j,i)==0
                obs.GPSL1C(j,i)=0;obs.GPSL2W(j,i)=0;obs.GPSC1W(j,i)=0;obs.GPSC2W(j,i)=0;
                EA.GPSel(j,i)=0;EA.GPSaz(j,i)=0;
                continue;
            end
            [el,az]=Get_EA(sx,sy,sz,gpsx(j,i)*1000,gpsy(j,i)*1000,gpsz(j,i)*1000);
            EA.GPSel(j,i)=el;
            EA.GPSaz(j,i)=az;
            if el<lim
                obs.GPSC1W(j,i)=0;obs.GPSC2W(j,i)=0;obs.GPSL1C(j,i)=0;obs.GPSL2W(j,i)=0;
                continue;
            end
        end
    end
end
% cut glonass obs
if ~isnan(find(strcmp(fields, 'GLOC1P' )))
    gloline=size(obs.GLOC1P,1);
    glosate=size(sate_mark.glo,2);
    if gloline<2880
        obs.GLOL1P(gloline+1:2880,:)=0;obs.GLOL2P(gloline+1:2880,:)=0;
        obs.GLOC1P(gloline+1:2880,:)=0;obs.GLOC2P(gloline+1:2880,:)=0;
    end
    glol=size(obs.GLOC1P,2);
    if glol<glosate
        obs.GLOL1P(:,glol+1:glosate)=0;obs.GLOL2P(:,glol+1:glosate)=0;
        obs.GLOC1P(:,glol+1:glosate)=0;obs.GLOC2P(:,glol+1:glosate)=0;
    end
    
    glodelete=find(sate_mark.glo==0);
    glonum=size(sate.glox,2);    
    if glonum<size(obs.GLOC1P,2)
        if ~isnan(glodelete)
            k=length(glodelete);
            for i=k:-1:1
                obs.GLOC1P(:,glodelete(i))=0;
                obs.GLOC2P(:,glodelete(i))=0;
                obs.GLOL1P(:,glodelete(i))=0;
                obs.GLOL2P(:,glodelete(i))=0;
            end
        end
        obs.GLOC1P=obs.GLOC1P(:,1:glonum);
        obs.GLOC2P=obs.GLOC2P(:,1:glonum);
        obs.GLOL1P=obs.GLOL1P(:,1:glonum);
        obs.GLOL2P=obs.GLOL2P(:,1:glonum);
    end
    obs.GLOC1P(isnan(obs.GLOC1P))=0;
    obs.GLOC2P(isnan(obs.GLOC2P))=0;
    obs.GLOL1P(isnan(obs.GLOL1P))=0;
    obs.GLOL2P(isnan(obs.GLOL2P))=0;
    
    glox=sate.glox;gloy=sate.gloy;gloz=sate.gloz;

    for i=1:glonum
        for j=1:2880
            if obs.GLOL1P(j,i)==0 || obs.GLOL2P(j,i)==0 || obs.GLOC1P(j,i)==0 || obs.GLOC2P(j,i)==0
                obs.GLOL1P(j,i)=0;obs.GLOL2P(j,i)=0;obs.GLOC1P(j,i)=0;obs.GLOC2P(j,i)=0;
                EA.GLOel(j,i)=0;EA.GLOaz(j,i)=0;
                continue;
            end
            [el,az]=Get_EA(sx,sy,sz,glox(j,i)*1000,gloy(j,i)*1000,gloz(j,i)*1000);
            EA.GLOel(j,i)=el;EA.GLOaz(j,i)=az;
            if el<lim
                obs.GLOC1P(j,i)=0;obs.GLOC2P(j,i)=0;obs.GLOL1P(j,i)=0;obs.GLOL2P(j,i)=0;
                continue;
            end
        end
    end
end
% cut bds obs
if ~isnan(find(strcmp(fields, 'BDSC2I' )))
    bdsline=size(obs.BDSC2I,1);
    bdssate=size(sate_mark.bds,2);
    if bdsline<2880
        obs.BDSL2I(bdsline+1:2880,:)=0;obs.BDSL6I(bdsline+1:2880,:)=0;
        obs.BDSC2I(bdsline+1:2880,:)=0;obs.BDSC6I(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSC2I,2);
    if bdsl<bdssate
        obs.BDSL2I(:,bdsl+1:bdssate)=0;obs.BDSL6I(:,bdsl+1:bdssate)=0;
        obs.BDSC2I(:,bdsl+1:bdssate)=0;obs.BDSC6I(:,bdsl+1:bdssate)=0;
    end
    
    bdsdelete=find(sate_mark.bds==0);
    bdsnum=size(sate.bdsx,2);
    if bdsnum<size(obs.BDSC2I,2)
        if ~isnan(bdsdelete)
            k=length(bdsdelete);
            for i=k:-1:1
                obs.BDSC2I(:,bdsdelete(i))=0;
                obs.BDSC6I(:,bdsdelete(i))=0;
                obs.BDSL2I(:,bdsdelete(i))=0;
                obs.BDSL6I(:,bdsdelete(i))=0;
            end
        end
        obs.BDSC2I=obs.BDSC2I(:,1:bdsnum);
        obs.BDSC6I=obs.BDSC6I(:,1:bdsnum);
        obs.BDSL2I=obs.BDSL2I(:,1:bdsnum);
        obs.BDSL6I=obs.BDSL6I(:,1:bdsnum);
    end
    obs.BDSC2I(isnan(obs.BDSC2I))=0;
    obs.BDSC6I(isnan(obs.BDSC6I))=0;
    obs.BDSL2I(isnan(obs.BDSL2I))=0;
    obs.BDSL6I(isnan(obs.BDSL6I))=0;
    
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;

    for i=1:bdsnum
        for j=1:2880
            if obs.BDSL2I(j,i)==0 || obs.BDSL6I(j,i)==0 || obs.BDSC2I(j,i)==0 || obs.BDSC6I(j,i)==0
                obs.BDSL2I(j,i)=0;obs.BDSL6I(j,i)=0;obs.BDSC2I(j,i)=0;obs.BDSC6I(j,i)=0;
                EA.BDSel(j,i)=0;EA.BDSaz(j,i)=0;
                continue;
            end
            [el,az]=Get_EA(sx,sy,sz,bdsx(j,i)*1000,bdsy(j,i)*1000,bdsz(j,i)*1000);
            EA.BDSel(j,i)=el;EA.BDSaz(j,i)=az;
            if el<lim
                obs.BDSC2I(j,i)=0;obs.BDSC6I(j,i)=0;obs.BDSL2I(j,i)=0;obs.BDSL6I(j,i)=0;
                continue;
            end
        end
    end
elseif  ~isnan(find(strcmp(fields, 'BDSC1I' )))
    bdsline=size(obs.BDSC1I,1);
    bdssate=size(sate_mark.bds,2);
    if bdsline<2880
        obs.BDSL1I(bdsline+1:2880,:)=0;obs.BDSL7I(bdsline+1:2880,:)=0;
        obs.BDSC1I(bdsline+1:2880,:)=0;obs.BDSC7I(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSC1I,2);
    if bdsl<bdssate
        obs.BDSL1I(:,bdsl+1:bdssate)=0;obs.BDSL7I(:,bdsl+1:bdssate)=0;
        obs.BDSC1I(:,bdsl+1:bdssate)=0;obs.BDSC7I(:,bdsl+1:bdssate)=0;
    end
    
    bdsdelete=find(sate_mark.bds==0);
    bdsnum=size(sate.bdsx,2);
    if bdsnum<size(obs.BDSC1I,2)
        if ~isnan(bdsdelete)
            k=length(bdsdelete);
            for i=k:-1:1
                obs.BDSC1I(:,bdsdelete(i))=0;
                obs.BDSC7I(:,bdsdelete(i))=0;
                obs.BDSL1I(:,bdsdelete(i))=0;
                obs.BDSL7I(:,bdsdelete(i))=0;
            end
        end
        obs.BDSC1I=obs.BDSC1I(:,1:bdsnum);
        obs.BDSC7I=obs.BDSC7I(:,1:bdsnum);
        obs.BDSL1I=obs.BDSL1I(:,1:bdsnum);
        obs.BDSL7I=obs.BDSL7I(:,1:bdsnum);
    end
    obs.BDSC1I(isnan(obs.BDSC1I))=0;
    obs.BDSC7I(isnan(obs.BDSC7I))=0;
    obs.BDSL1I(isnan(obs.BDSL1I))=0;
    obs.BDSL7I(isnan(obs.BDSL7I))=0;
    
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;

    for i=1:bdsnum
        for j=1:2880
            if obs.BDSL1I(j,i)==0 || obs.BDSL7I(j,i)==0 || obs.BDSC1I(j,i)==0 || obs.BDSC7I(j,i)==0
                obs.BDSL1I(j,i)=0;obs.BDSL7I(j,i)=0;obs.BDSC1I(j,i)=0;obs.BDSC7I(j,i)=0;
                EA.BDSel(j,i)=0;EA.BDSaz(j,i)=0;
                continue;
            end
            [el,az]=Get_EA(sx,sy,sz,bdsx(j,i)*1000,bdsy(j,i)*1000,bdsz(j,i)*1000);
            EA.BDSel(j,i)=el;EA.BDSaz(j,i)=az;
            if el<lim
                obs.BDSC1I(j,i)=0;obs.BDSC7I(j,i)=0;obs.BDSL1I(j,i)=0;obs.BDSL7I(j,i)=0;
                continue;
            end
        end
    end
end
% cut galileo obs
if ~isnan(find(strcmp(fields, 'GALC1X' )))
    galxline=size(obs.GALC1X,1);
    galsate=size(sate_mark.gal,2);
    
    if galxline<2880
        obs.GALL1X(galxline+1:2880,:)=0;obs.GALL5X(galxline+1:2880,:)=0;
        obs.GALC1X(galxline+1:2880,:)=0;obs.GALC5X(galxline+1:2880,:)=0;
    end
    gall=size(obs.GALC1X,2);
    if gall<galsate
        obs.GALL1X(:,gall+1:galsate)=0;obs.GALL5X(:,gall+1:galsate)=0;
        obs.GALC1X(:,gall+1:galsate)=0;obs.GALC5X(:,gall+1:galsate)=0;
    end
    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);
    if galnum<size(obs.GALC1X,2)
        if ~isnan(galdelete)
            k=length(galdelete);
            for i=k:-1:1
                obs.GALC1X(:,galdelete(i))=0;
                obs.GALC5X(:,galdelete(i))=0;
                obs.GALL1X(:,galdelete(i))=0;
                obs.GALL5X(:,galdelete(i))=0;
            end
        end
        obs.GALC1X=obs.GALC1X(:,1:galnum);
        obs.GALC5X=obs.GALC5X(:,1:galnum);
        obs.GALL1X=obs.GALL1X(:,1:galnum);
        obs.GALL5X=obs.GALL5X(:,1:galnum);
    end
    obs.GALC1X(isnan(obs.GALC1X))=0;
    obs.GALC5X(isnan(obs.GALC5X))=0;
    obs.GALL1X(isnan(obs.GALL1X))=0;
    obs.GALL5X(isnan(obs.GALL5X))=0;
    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;

    for i=1:galnum
        for j=1:2880
            if obs.GALL1X(j,i)==0 || obs.GALL5X(j,i)==0 || obs.GALC1X(j,i)==0 || obs.GALC5X(j,i)==0
                obs.GALL1X(j,i)=0;obs.GALL5X(j,i)=0;obs.GALC1X(j,i)=0;obs.GALC5X(j,i)=0;
                EA.GALel(j,i)=0;EA.GALaz(j,i)=0;
                continue;
            end
            [el,az]=Get_EA(sx,sy,sz,galx(j,i)*1000,galy(j,i)*1000,galz(j,i)*1000);
            EA.GALel(j,i)=el;EA.GALaz(j,i)=az;
            if el<lim
                obs.GALC1X(j,i)=0;obs.GALC5X(j,i)=0;
                obs.GALL1X(j,i)=0;obs.GALL5X(j,i)=0;
                continue;
            end
        end
    end
elseif ~isnan(find(strcmp(fields, 'GALC1C' )))
    galline=size(obs.GALC1C,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALL1C(galline+1:2880,:)=0;obs.GALL5Q(galline+1:2880,:)=0;
        obs.GALC1C(galline+1:2880,:)=0;obs.GALC5Q(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC1C,2);
    if gall<galsate
        obs.GALL1C(:,gall+1:galsate)=0;obs.GALL5Q(:,gall+1:galsate)=0;
        obs.GALC1C(:,gall+1:galsate)=0;obs.GALC5Q(:,gall+1:galsate)=0;
    end
    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);
        
    if galnum<size(obs.GALC1C,2)
        if ~isnan(galdelete)
            k=length(galdelete);
            for i=k:-1:1
                obs.GALC1C(:,galdelete(i))=0;
                obs.GALC5Q(:,galdelete(i))=0;
                obs.GALL1C(:,galdelete(i))=0;
                obs.GALL5Q(:,galdelete(i))=0;
            end
        end
        obs.GALC1C=obs.GALC1C(:,1:galnum);
        obs.GALC5Q=obs.GALC5Q(:,1:galnum);
        obs.GALL1C=obs.GALL1C(:,1:galnum);
        obs.GALL5Q=obs.GALL5Q(:,1:galnum);
    end
    obs.GALC1C(isnan(obs.GALC1C))=0;
    obs.GALC5Q(isnan(obs.GALC5Q))=0;
    obs.GALL1C(isnan(obs.GALL1C))=0;
    obs.GALL5Q(isnan(obs.GALL5Q))=0;
    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;

    for i=1:galnum
        for j=1:2880
            if obs.GALL1C(j,i)==0 || obs.GALL5Q(j,i)==0 || obs.GALC1C(j,i)==0 || obs.GALC5Q(j,i)==0
                obs.GALL1C(j,i)=0;obs.GALL5Q(j,i)=0;
                obs.GALC1C(j,i)=0;obs.GALC5Q(j,i)=0;
                EA.GALel(j,i)=0;EA.GALaz(j,i)=0;
                continue;
            end
            [el,az]=Get_EA(sx,sy,sz,galx(j,i)*1000,galy(j,i)*1000,galz(j,i)*1000);
            EA.GALel(j,i)=el;EA.GALaz(j,i)=az;
            if el<lim  
                obs.GALC1C(j,i)=0;obs.GALC5Q(j,i)=0;
                obs.GALL1C(j,i)=0;obs.GALL5Q(j,i)=0;
                continue;
            end
        end
    end
end
end
%% 
function [E,A]= Get_EA(sx,sy,sz,x,y,z)
%GET_EL Summary of this function goes here
%   Detailed explanation goes here
if x==0 && y==0 && z==0
E=0;
A=0;
else
[sb,sl]=XYZtoBLH(sx,sy,sz); %station latitude and longitude
T=[-sin(sb)*cos(sl) -sin(sb)*sin(sl) cos(sb);
    -sin(sl)               cos(sl)         0;
    cos(sb)*cos(sl) cos(sb)*sin(sl)  sin(sb)];%transition matrix(XYZ to NEU)
deta_xyz=[x,y,z]-[sx,sy,sz];
NEU=T*(deta_xyz)';
E=atan(NEU(3)/sqrt(NEU(1)*NEU(1)+NEU(2)*NEU(2)));
A=atan(abs(NEU(2)/NEU(1)));
if NEU(1)>0
    if NEU(2)>0
    else
        A=2*pi-A;
    end
else
    if NEU(2)>0
        A=pi-A;
    else
        A=pi+A;
    end 
end
end
end
%% 
function index=find_sp3(list,doy)
%find related sp3
len=length(list);
doys=linspace(0,0,len);
for i=1:len
    doys(i)=str2double(list(i).name(3:7));
end
index=find(doys==str2double(doy),1);
end 


