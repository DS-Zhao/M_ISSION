function get_RegionalROTI(rotipath,lonlat)
%% Regional ROTI data calculation
% INPUT:
%     rotipath: This software calculates and saves the ROTI.mat data folder path
%     range: Grid division setting
%            Includes: minlat, maxlat, minlon, maxlon, dlat, dlon
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
lat_min = lonlat(1);
lat_max = lonlat(2);
lon_min = lonlat(3);
lon_max = lonlat(4);
lat_resolution = lonlat(5); %Latitude resolution, unit: degrees
lon_resolution = lonlat(6); %Longitude resolution, unit: degrees
num_lats = ceil((lat_max - lat_min) / lat_resolution);
num_lons = ceil((lon_max - lon_min) / lon_resolution);
%Initialize a cell array to store the grid for each time period
ROTIgrid = cell(1, 24);
list_roti = dir([rotipath, '/*.mat']);
grandPath=fileparts(fileparts(rotipath));
for i = 1:length(list_roti)
    rotiname = list_roti(i).name;
    doy=rotiname(5:9);
    if strcmp(rotiname(end-7:end-4), 'ROTI')
        load([rotipath, '\', rotiname]);
        obsPath = fullfile(grandPath, 'raw_OBS_cut', rotiname(5:9), [rotiname(1:9), '.mat']);
        if exist(obsPath, 'file')
            load(obsPath);

            fields2 = fieldnames(EA);
            sx = coor(1); sy = coor(2); sz = coor(3);

            if any(strcmp(fields2, 'GPSel')) && strcmp(rotiname(end-10:end-8), 'GPS')
                [~, GPSB, GPSL] = Get_dif(sx, sy, sz, EA.GPSel, EA.GPSaz);
                L = mod(GPSL * 180 / pi + 180, 360) - 180;
                B = GPSB * 180 / pi;
                ROTI=GPSROTI;

            elseif any(strcmp(fields2, 'GLOel')) && strcmp(rotiname(end-10:end-8), 'GLO')
                [~, GLOB, GLOL] = Get_dif(sx, sy, sz, EA.GLOel, EA.GLOaz);
                L = mod(GLOL * 180 / pi + 180, 360) - 180;
                B = GLOB * 180 / pi;
                ROTI=GLOROTI;

            elseif any(strcmp(fields2, 'GALel')) && strcmp(rotiname(end-10:end-8), 'GAL')
                [~, GALB, GALL] = Get_dif(sx, sy, sz, EA.GALel, EA.GALaz);
                L = mod(GALL * 180 / pi + 180, 360) - 180;
                B = GALB * 180 / pi;
                ROTI=GALROTI;

            elseif any(strcmp(fields2, 'GALel')) && strcmp(rotiname(end-11:end-8), 'GALX')
                [~, GALB, GALL] = Get_dif(sx, sy, sz, EA.GALel, EA.GALaz);
                L = mod(GALL * 180 / pi + 180, 360) - 180;
                B = GALB * 180 / pi;
                ROTI=GALXROTI;

            elseif any(strcmp(fields2, 'BDSel')) && strcmp(rotiname(end-10:end-8), 'BDS')
                [~, BDSB, BDSL] = Get_dif(sx, sy, sz, EA.BDSel, EA.BDSaz);
                L = mod(BDSL * 180 / pi + 180, 360) - 180;
                B = BDSB * 180 / pi;
                ROTI=BDSROTI;
            end
            %Process data for each time period
            ROTIgrid=hourgrid(L,B,ROTI,num_lats,num_lons,lonlat,ROTIgrid);
        end
    end
end
count = sum(~cellfun(@isempty, ROTIgrid));
if count==24
    fPath=[fileparts(rotipath),'\regionalROTI\'];
    if exist(fPath,'dir')==0
        mkdir(fPath);
    end
    ffPath=[fPath ,doy rotiname(end-10:end-8) 'regionalROTI.mat'];
    save(ffPath,'ROTIgrid','lonlat','-mat');
end
end