function [] = hourROTI(LL, sub, minlon, minlat, maxlon, maxlat,ax2)
%% Draw a one hour regional map
%INPUT:
%     LL:Grid data of ROTI
%     sub：Current drawing hour
%     minlon,minlat,maxlon,maxlat：The latitude and longitude intervals drawn
%     ax2: Drawing window parameters
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
%Range of image data acquisition
data_min = min(LL{1,sub}(:), [], 'all', 'omitnan'); %Ignore NaN values using omitnan
data_max = max(LL{1,sub}(:), [], 'all', 'omitnan');

%Create geogrid data reference objects
R = georasterref('RasterSize', size(LL{1,sub}), 'Latlim', [minlat maxlat], 'Lonlim', [minlon maxlon]);

%Set the map projection method
ax = axesm('MapProjection','robinson','MapLatLimit',[minlat maxlat],'MapLonLimit',[minlon maxlon]);
axis off; 

%Set map grid and labels
setm(ax, 'GLineStyle', '--', 'Frame', 'on', 'GlineWidth', 1);
meridian_interval = round(abs(maxlon-minlon)/5);
parallel_interval = round(abs(maxlat-minlat)/5);
setm(ax, ...
    'MlabelLocation', meridian_interval, ...
    'PlabelLocation', parallel_interval, ...
    'MeridianLabel', 'on', ...
    'ParallelLabel', 'on', ...
    'MlineLocation', meridian_interval, ...
    'PlineLocation', parallel_interval, ...
    'MLabelParallel', 'south' ...
    );
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [250,227,158]./255);

%Display geographic data
h = geoshow(LL{1,sub}, R, 'DisplayType', 'texturemap');
h.AlphaDataMapping = 'none'; 
h.FaceAlpha = 'texturemap';
alpha(h, double(~isnan(LL{1,sub}))); 

%Set coordinate axis box
set(gca, 'box', 'on');
colormap(ax2,"parula") ;

%Set Title
title(['ROTI Hour - ', num2str(sub)]);

%Create a color bar and configure it as a gradient color
cb = colorbar;
cb.Location = 'eastoutside'; %Set the position of the color bar
cb.TickLength = 0.02; %Set the length of the color bar scale

%Set the range of the color bar based on the data range
clim([0, data_max]); %Set the range of color bars

end