function []=Get_DIXSG(grandPath,Sites_Info,lev,first,step,minlon,minlat,maxlon,maxlat,dlat,dlon,sys)
%% Calculate DIXSG data
% INPUT:
%     grandPath: main path
%     Sites_Info: name and coordinate information of the stations
%     lev:Number of Grades
%     first:first-degree
%     step:step lengt
%     minlon,minlat,maxlon,maxlat,dlat,dlon:Regional longitude and latitude range and grid division size
%     sys: Satellite system selection parameter status
%% written by Zhang P,Zhang X. et al., 2024/08
%% --------------------------------------------------------------
stations=Sites_Info.name;
doys=Sites_Info.doy;
current_doy=num2str(unique(doys));
cROTfoldername=dir(fullfile(grandPath,'ivcROT'));
subfolderNames = {};
MBL=[maxlon,minlon,maxlat,minlat,lev];
for i = 1:length(cROTfoldername)
    if cROTfoldername(i).isdir && ~strcmp(cROTfoldername(i).name, '.') && ~strcmp(cROTfoldername(i).name, '..')
        subfolderNames = [subfolderNames, cROTfoldername(i).name];
    end
end
cROTname2={};
for i = 1:length(subfolderNames)
    folderName = subfolderNames(i);
    folderName = char(folderName);
    folderNumber = folderName(end-4:end);
    if strcmp(folderNumber, current_doy)
        cROTname2{end+1} = folderName;
    end
end

for i = 1:numel(cROTname2)
    if strncmp(cROTname2{i}, 'GPSc', 4) && sys(1)~=0
        cROTpath=[grandPath '/ivcROT' '/GPScROT' current_doy];
        D_Resultname=dir([cROTpath,'/*.mat']);
        D_Resultname=now_files(stations,D_Resultname);
        num_files = numel(D_Resultname);
        if num_files < 2
            continue;
        end
        satellite='GPS';
        douSITE=Get_TALL(D_Resultname,cROTpath,lev,first,step,minlon,minlat,maxlon,maxlat,dlat,dlon,satellite,current_doy);
        DIXSG(grandPath,satellite,current_doy,douSITE,MBL);
    elseif strncmp(cROTname2{i}, 'GLOc', 4) && sys(2)~=0
        cROTpath=[grandPath '/ivcROT' '/GLOcROT' current_doy];
        D_Resultname=dir([cROTpath,'/*.mat']);
        D_Resultname=now_files(stations,D_Resultname);
        num_files = numel(D_Resultname);
        if num_files < 2
            continue;
        end
        satellite='GLO';
        douSITE=Get_TALL(D_Resultname,cROTpath,lev,first,step,minlon,minlat,maxlon,maxlat,dlat,dlon,satellite,current_doy);
        DIXSG(grandPath,satellite,current_doy,douSITE,MBL);
    elseif strncmp(cROTname2{i}, 'BDSc', 4) && sys(4)~=0
        cROTpath=[grandPath '/ivcROT' '/BDScROT' current_doy];
        D_Resultname=dir([cROTpath,'/*.mat']);
        D_Resultname=now_files(stations,D_Resultname);
        num_files = numel(D_Resultname);
        if num_files < 2
            continue;
        end
        satellite='BDS';
        douSITE=Get_TALL(D_Resultname,cROTpath,lev,first,step,minlon,minlat,maxlon,maxlat,dlat,dlon,satellite,current_doy);
        DIXSG(grandPath,satellite,current_doy,douSITE,MBL);
    elseif strncmp(cROTname2{i}, 'GALc', 4) && sys(3)~=0
        cROTpath=[grandPath '/ivcROT' '/GALcROT' current_doy];
        D_Resultname=dir([cROTpath,'/*.mat']);
        D_Resultname=now_files(stations,D_Resultname);
        num_files = numel(D_Resultname);
        if num_files < 2
            continue;
        end
        satellite='GAL';
        douSITE=Get_TALL(D_Resultname,cROTpath,lev,first,step,minlon,minlat,maxlon,maxlat,dlat,dlon,satellite,current_doy);
        DIXSG(grandPath,satellite,current_doy,douSITE,MBL);
    end
end
end
%%
function douSITE=Get_TALL(D_Resultname,cROTpath,lev,first,step,minlon,minlat,maxlon,maxlat,dlat,dlon,satellite,current_doy)
grandPath=fileparts(fileparts(cROTpath));
len=length(D_Resultname);
douSITE={};
for msta=1:len-1
    for sta=msta+1:len
        load([cROTpath,'/',D_Resultname(msta).name],'-mat');
        if strcmp(satellite, 'GPS')
            F1=GPS_B_L_cROT;
        elseif strcmp(satellite, 'GLO')
            F1=GLO_B_L_cROT;
        elseif strcmp(satellite, 'BDS')
            F1=BDS_B_L_cROT;
        elseif strcmp(satellite, 'GAL')
            F1=GAL_B_L_cROT;
        elseif strcmp(satellite, 'GALX')
            F1=GALX_B_L_cROT;
        end
        SITE{1}=upper(D_Resultname(msta).name(1:4));
        load([cROTpath,'/',D_Resultname(sta).name],'-mat');
        if strcmp(satellite, 'GPS')
            F2=GPS_B_L_cROT;
        elseif strcmp(satellite, 'GLO')
            F2=GLO_B_L_cROT;
        elseif strcmp(satellite, 'BDS')
            F2=BDS_B_L_cROT;
        elseif strcmp(satellite, 'GAL')
            F2=GAL_B_L_cROT;
        elseif strcmp(satellite, 'GALX')
            F2=GALX_B_L_cROT;
        end
        SITE{2}=upper(D_Resultname(sta).name(1:4));
        %---------------d-----------------------
        fname=[grandPath '\' 'ivTall\' satellite current_doy '\' SITE{1} '-' SITE{2} '.mat'];
        % if isfile(fname)
        %     douSITE{end+1}=[SITE{1} '-' SITE{2}];
        %     fclose all;
        %     continue;
        % end
        D=1000; %km
        g=min(size(F1.B,2),size(F2.B,2));
        d=cal_DIPP2(F1,F2);
        if any(d(:)>1000)||all(isnan(d(:)))
            continue;
        end
        level = zeros(1, lev);
        for i = 1:lev
            level(i) = first + (i - 1) * step;
        end
        Level = cell(1, numel(level));
        for i = 1:numel(level)
            Level{i} = ['d' num2str(level(i))];
        end

        FDIXSG=0;
        for i=1:lev
            [F1.cROT, F2.cROT] = resizeArraysToFit(F1.cROT, F2.cROT);
            DIXSG.(Level{1,i})=(abs(F1.cROT-F2.cROT)/level(i)).^3.*(d/D).^-1;
            DIXSG.(Level{1,i})(DIXSG.(Level{1,i})>1)=1;
            DIXSG.(Level{1,i})(DIXSG.(Level{1,i})<1)=0;
            FDIXSG=DIXSG.(Level{1,i})+FDIXSG; % Simplify each sensitivity data to>1=1,<1=0
        end
        if all(isnan(FDIXSG(:)))
            continue;
        end
        %---------------getCenterLonLat-----------------------
        F1.L=mod(F1.L* 180 / pi + 180, 360) - 180;
        F2.L=mod(F2.L* 180 / pi + 180, 360) - 180; %Turn to angle
        F1.B=F1.B*180/pi  ; F2.B=F2.B*180/pi; %Turn to angle
        IPP_pair.lon=zeros(2880,g);
        IPP_pair.lat=zeros(2880,g);
        for i=1:g
            for j=1:2880
                if F1.L(j,i)~=0 && F2.L(j,i)~=0 && F1.B(j,i)~=0 && F2.B(j,i)~=0
                    lon1 = F1.L(j, i);lat1 = F1.B(j, i);lon2 = F2.L(j, i);lat2 = F2.B(j, i);
                    blon = (lon1 + lon2) / 2;
                    blat = (lat1 + lat2) / 2;
                end
                IPP_pair.lon(j,i)=blon;
                IPP_pair.lat(j,i)=blat;
            end
        end

        %---------------Partition grid-----------------------

        numlon = ceil((maxlon-minlon)/dlon);%Longitude grid
        numlat = ceil((maxlat-minlat)/dlat);%Latitude grid
        H = cell(1, numlat);
        for i = 1:numlat
            H{i} = ['L' num2str(i)];
        end
        Grids=Get_grid(IPP_pair,FDIXSG,minlon,minlat,maxlon,maxlat,numlon,numlat,dlon,dlat,H,g);
        if isempty(Grids)
            continue;
        end
        Tall=cell(1,24);
        Tmax=[];
        for k=1:24
            for l=1:numlon
                Max=NaN(numlat,1);
                data=Grids{1,k}{1,l};

                if ~isempty(data)
                    for h=numlat:-1:1
                        if isfield(data,H{1,h})
                            Max(h,1)=max(data.(H{1,h}));
                        end
                    end
                end
                Tmax=[Tmax,Max]; %Total grid data within 1 hour
            end
            Tall{1,k}=Tmax; %Total grid data within 24 hours 1×24
            Tmax=[];
        end
        if exist([grandPath '\' 'ivTall\' satellite current_doy '\'  ],'dir')==0
            mkdir([grandPath '\' 'ivTall\' satellite current_doy '\' ]);
        end
        fname=[grandPath '\' 'ivTall\' satellite current_doy '\' SITE{1} '-' SITE{2} '.mat'];
        douSITE{end+1}=[SITE{1} '-' SITE{2}];
        save(fname,'Tall','-mat');
        fclose all;
    end
end
end
function DIXSG(grandPath,satellite,current_doy,list_tall,MBL)
LL=cell(1,24);
DIXSG=zeros(1,24);
tic
for k=1:24
    com=0;
    Num=0;
    tran=0;
    for i=1:size(list_tall,2) %At least two sets of data
        A=load([grandPath '\' 'ivTall\' satellite current_doy '\' list_tall{i},'.mat']);
        a=A.Tall{1,k};
        if ~all(all(isnan(a))) %If all observations in a pair are NaN, discard them
            [m,n]=size(a);
            num=m *n-numel(find(isnan(a)));
            Num=Num+num;
            com=com+sum(sum(a,'omitnan'));

            if tran==0
                tran=a;
            else
                for I=1:m
                    for J=1:n
                        tran(I,J)=max(tran(I,J),a(I,J));
                    end
                end
            end
        end
    end
    aDIXSG(k)=com/Num;
    LL{1,k}=tran; %DIXS every hour
end
if exist([grandPath '\' 'resDIXSG\' satellite 'DIXSG' current_doy '\' ],'dir')==0
    mkdir([grandPath '\' 'resDIXSG\' satellite 'DIXSG' current_doy ]);
end
fname=[grandPath '\' 'resDIXSG\' satellite 'DIXSG' current_doy '\' satellite current_doy 'DIXSG' '.mat'];
save(fname,'aDIXSG','LL','MBL','-mat');
fclose all;
end
%%Observation file corresponding file
function new_files=now_files(stations,files)
file_features = cell(length(files)-2, 1);
for i = 3:length(files)
    file_name = files(i).name;
    file_features{i-2} = file_name(1:4);  %Extract the first four feature names
end
%Check if each feature name exists in stations
valid_files_indices = [];
for i = 1:length(file_features)
    if ismember(file_features{i}, stations)
        valid_files_indices = [valid_files_indices, i+2];
    end
end
new_files = files(valid_files_indices);
end
%% Uniform array size
function [A_resized, B_resized] = resizeArraysToFit(A, B)

[rowsA, colsA] = size(A);
[rowsB, colsB] = size(B);

minRows = min(rowsA, rowsB);
minCols = min(colsA, colsB);

A_resized = A(1:minRows, 1:minCols);
B_resized = B(1:minRows, 1:minCols);
end