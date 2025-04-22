function [] = read_sp3( s_ipath,s_opath,sate_mark)
%% read sp3 files to struct
% INPUT:
%     s_ipath: storage path of *.sp3 files
%     s_opath: storage path of satellite coordinate files
%     sate_mark: satellite status identification
% SAVE:
%      */SP3/doysp3.mat:Sp3 data
%% written by Jin R et al., 2012/5/20, doi:10.1007/s10291-012-0279-3
%% modified by Zhou C. et al., 2021/12/12
%% modified by Zhang P. et al., 2024/08
%% --------------------------------------------------------------------
list_obs=dir([s_ipath '/*.sp3']);
len=length(list_obs);
if len<3
 error('need at least three SP3 files !!!');
end
for i=1:len-2
    GN=list_obs(i+1).name(12:18);
    filename=strcat(GN,'sp3.mat');
    % if isfile([s_opath,'/',filename])
    %     continue;
    % else   
    pr_obs=[s_ipath '/' list_obs(i).name];  
    cu_obs=[s_ipath '/' list_obs(i+1).name];
    nx_obs=[s_ipath '/' list_obs(i+2).name];
     [gpsx1,gpsy1,gpsz1,glox1,gloy1,gloz1,galx1,galy1,galz1,bdsx1,bdsy1,bdsz1]=r_sp3(pr_obs);
    [gpsx2,gpsy2,gpsz2,glox2,gloy2,gloz2,galx2,galy2,galz2,bdsx2,bdsy2,bdsz2]=r_sp3(cu_obs);
    [gpsx3,gpsy3,gpsz3,glox3,gloy3,gloz3,galx3,galy3,galz3,bdsx3,bdsy3,bdsz3]=r_sp3(nx_obs);
    numgps=min([size(gpsx1,2),size(gpsx2,2),size(gpsx3,2)]);
    numglo=min([size(glox1,2),size(glox2,2),size(glox3,2)]);
    numgal=min([size(galx1,2),size(galx2,2),size(galx3,2)]);
    numbds=min([size(bdsx1,2),size(bdsx2,2),size(bdsx3,2)]);
    sate.gpsx=[];sate.gpsy=[];sate.gpsz=[];
    sate.glox=[];sate.gloy=[];sate.gloz=[];
    sate.galx=[];sate.galy=[];sate.galz=[];
    sate.bdsx=[];sate.bdsy=[];sate.bdsz=[];
    if numgps~=0
    [sate.gpsx,sate.gpsy,sate.gpsz]=interplotation(numgps,gpsx1,gpsy1,gpsz1,gpsx2,gpsy2,gpsz2,gpsx3,gpsy3,gpsz3);
    lsate = size(sate.gpsx, 2);
    if lsate < size(sate_mark.gps, 2)
        sate_mark2 = sate_mark.gps(1, 1:lsate);
    else
        sate_mark2 = sate_mark.gps;
    end
    gpsdelete=find(sate_mark2==0);
    if ~isnan(gpsdelete)
        kg=length(gpsdelete);
        for j=kg:-1:1
            sate.gpsx(:,gpsdelete(j))=0;
            sate.gpsy(:,gpsdelete(j))=0;
            sate.gpsz(:,gpsdelete(j))=0;
        end
    end
    end
    if numglo~=0
    [sate.glox,sate.gloy,sate.gloz]=interplotation(numglo,glox1,gloy1,gloz1,glox2,gloy2,gloz2,glox3,gloy3,gloz3);
    lsate = size(sate.glox, 2);
    if lsate < size(sate_mark.glo, 2)
        sate_mark2 = sate_mark.glo(1, 1:lsate);
    else
        sate_mark2 = sate_mark.glo;
    end
    glodelete=find(sate_mark2==0);
    if ~isnan(glodelete)
        kr=length(glodelete);
        for j=kr:-1:1
            sate.glox(:,glodelete(j))=0;
            sate.gloy(:,glodelete(j))=0;
            sate.gloz(:,glodelete(j))=0;
        end
    end
    end
    if numgal~=0
    [sate.galx,sate.galy,sate.galz]=interplotation(numgal,galx1,galy1,galz1,galx2,galy2,galz2,galx3,galy3,galz3);
    lsate = size(sate.galx, 2);
    if lsate < size(sate_mark.gal, 2)
        sate_mark2 = sate_mark.gal(1, 1:lsate);
    else
        sate_mark2 = sate_mark.gal;
    end
    galdelete=find(sate_mark2==0);
    if ~isnan(galdelete)
        ke=length(galdelete);
        for j=ke:-1:1
            sate.galx(:,galdelete(j))=0;
            sate.galy(:,galdelete(j))=0;
            sate.galz(:,galdelete(j))=0;
        end
    end
    end
    if numbds~=0
    [sate.bdsx,sate.bdsy,sate.bdsz]=interplotation(numbds,bdsx1,bdsy1,bdsz1,bdsx2,bdsy2,bdsz2,bdsx3,bdsy3,bdsz3);
    lsate = size(sate.bdsx, 2);
    if lsate < size(sate_mark.bds, 2)
        sate_mark2 = sate_mark.bds(1, 1:lsate);
    else
        sate_mark2 = sate_mark.bds;
    end

    bdsdelete = find(sate_mark2 == 0);
    if ~isempty(bdsdelete)
        kb = length(bdsdelete);
        for j = kb:-1:1
            sate.bdsx(:, bdsdelete(j)) = 0;
            sate.bdsy(:, bdsdelete(j)) = 0;
            sate.bdsz(:, bdsdelete(j)) = 0;
        end
    end
    end
    % GN=list_obs(i+1).name(1:7);
    % filename=strcat(GN,'sp3.mat');
    if exist(s_opath,'dir')==0
        mkdir(s_opath);
    end
    if ~isempty(sate.gpsx)
    save([s_opath,'/',filename],'sate','-mat');
    end
    % end
end
end
%% ----------------subfunction-----------------
function [GPSX,GPSY,GPSZ,GLOX,GLOY,GLOZ,GALX,GALY,GALZ,BDSX,BDSY,BDSZ] = r_sp3(path)
GPSX=[];GPSY=[];GPSZ=[];GLOX=[];GLOY=[];GLOZ=[];GALX=[];GALY=[];GALZ=[];BDSX=[];BDSY=[];BDSZ=[];
fid = fopen(path,'r');
line=fgetl(fid);
ep=0;%epoch number
day=line(12:13);
hm=5;
while 1
    if ~ischar(line), break, end
    line=fgetl(fid);
    if strcmp(line(1),'#')&& strcmp(line(27:29),'900')
     hm=15;
    end
    if strcmp(line(1),'*') && strcmp(line(12:13),day)
        h=str2double(line(15:16));
        m=str2double(line(18:19));
        ep=h*(60/hm)+round(m/hm)+1;
        continue;
    end
    if strcmp(line(1),'*') && ~strcmp(line(12:13),day)
        break;
    end
    if length(line)>1 && strcmp(line(1:2),'PG')
        sv=str2double(line(3:4));
        GPSX(ep,sv) = str2double(line(5:18));
        GPSY(ep,sv) = str2double(line(19:32));
        GPSZ(ep,sv) = str2double(line(33:46));
        continue;
    end
    if length(line)>1 && strcmp(line(1:2),'PR')
        sv=str2double(line(3:4));
        GLOX(ep,sv) = str2double(line(5:18));
        GLOY(ep,sv) = str2double(line(19:32));
        GLOZ(ep,sv) = str2double(line(33:46));
        continue;
    end
    if length(line)>1 && strcmp(line(1:2),'PE')
        sv=str2double(line(3:4));
        GALX(ep,sv) = str2double(line(5:18));
        GALY(ep,sv) = str2double(line(19:32));
        GALZ(ep,sv) = str2double(line(33:46));
        continue;
    end
    if length(line)>1 && strcmp(line(1:2),'PC')
        sv=str2double(line(3:4));
        BDSX(ep,sv) = str2double(line(5:18));
        BDSY(ep,sv) = str2double(line(19:32));
        BDSZ(ep,sv) = str2double(line(33:46));
        continue;
    end
    
end
fclose(fid);
end
%% ----------------subfunction-------------------
function [interp_x2,interp_y2,interp_z2]=interplotation(satenum,x1,y1,z1,x2,y2,z2,x3,y3,z3)
interp_x2=zeros(2880,satenum);
interp_y2=zeros(2880,satenum);
interp_z2=zeros(2880,satenum);
[x1, x2, x3] = align_dimensions(x1, x2, x3);
[y1, y2, y3] = align_dimensions(y1, y2, y3);
[z1, z2, z3] = align_dimensions(z1, z2, z3);
if length(x2)>96
x2=[x1(end-3:end,:);x2;x3(1:5,:)];
y2=[y1(end-3:end,:);y2;y3(1:5,:)];
z2=[z1(end-3:end,:);z2;z3(1:5,:)];
dpp=9;
dt1=4;dt2=5;dp=10;
else
x2=[x1(end-1:end,:);x2;x3(1:3,:)];
y2=[y1(end-1:end,:);y2;y3(1:3,:)];
z2=[z1(end-1:end,:);z2;z3(1:3,:)];
dt1=2;dt2=3;dp=30;dpp=5;
end 
NUMCOL=length(x2)-dpp;
% m_t=linspace(-40,2880+40,length(x2));
m_t=linspace(-dt1*dp,2880+dt1*dp,length(x2));
for i=1:satenum
    for j=1:NUMCOL     
        tt=m_t(j:j+dpp);
        x=x2((j:dpp+j),i)';y=y2((j:dpp+j),i)';z=z2((j:dpp+j),i)';
        t0=linspace(m_t(j+dt1),m_t(j+dt2)-1,dp);
        interp_x2((dp*j-(dp-1)):dp*j,i)=interp_lag(tt,x,t0)';
        interp_y2((dp*j-(dp-1)):dp*j,i)=interp_lag(tt,y,t0)';
        interp_z2((dp*j-(dp-1)):dp*j,i)=interp_lag(tt,z,t0)';
    end
end
end
%% ----------------subfunction----------------
function y0 = interp_lag (x, y, x0)
n=length(x);
y0=zeros(size(x0));
for k=1:n
    t=1;
    for i=1:n
        if i~=k
            t=t.*(x0-x(i))/(x(k)-x(i));
        end
    end
    y0=y0+t*y(k);
end
end
%% ----------------subfunction-----------------
function [s1_out, s2_out, s3_out] = align_dimensions(s1, s2, s3)
    % Input:
    %   s1, s2, s3: Three matrices input
    % Output:
    %   s1_out, s2_out, s3_out: Output data after dimension adjustment
%% written by Zhang P. et al., 2024/08

    % Get the number of columns for each input matrix (assuming dimension alignment of columns)
    [ds1,dims_s1] = size(s1);
    [ds2,dims_s2] = size(s2);
    [ds3,dims_s3] = size(s3);
    
    if any([ds1,ds2,ds3]>96) && any([ds1,ds2,ds3]==96)
        
            if ds1>96
                s1=s1(1:3:288,:);
            end
            if ds2>96
                s2=s2(1:3:288,:);
            end
            if ds3>96
                s3=s3(1:3:288,:);
            end
    end

    % Find the minimum number of columns
    min_cols = min([dims_s1, dims_s2, dims_s3]);

    % Adjust the number of columns in each matrix to the minimum number of columns
    s1_out = s1(:, 1:min_cols);
    s2_out = s2(:, 1:min_cols);
    s3_out = s3(:, 1:min_cols);
end
