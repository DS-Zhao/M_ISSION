function []=hourDIXSG(lev,LL,sub,minlon,minlat,maxlon,maxlat,ax1)
%% Draw a one hour regional map
%INPUT:
%     lev: sensitivity
%     LL:Grid data of DIXSG
%     sub：Current drawing hour
%     minlon,minlat,maxlon,maxlat：The latitude and longitude intervals drawn
%     ax1: Drawing window parameters
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
R = georasterref('RasterSize', size(LL{1,sub}), 'Latlim', [minlat maxlat], 'Lonlim', [minlon maxlon]);%Geographic Grid Data Reference Object (Class)
ax = axesm('MapProjection','robinson','MapLatLimit',[minlat maxlat],'MapLonLimit',[minlon maxlon]);%Set the equidistant cylindrical projection method for the map
axis off%Close the local coordinate axis system
setm(ax,'GLineStyle','--' ,'Frame','on','GlineWidth',1)%Specify grid line shape, draw frame 
setm(ax,...
    'MlabelLocation',round(abs(maxlon-minlon)/5),...
    'PlabelLocation',round(abs(maxlat-minlat)/5),...%Draw latitude scale labels only at specified values
    'MeridianLabel','on',...%Display longitude scale label
    'ParallelLabel','on',...%Display latitude scale labels
    'MlineLocation',round(abs(maxlon-minlon)/5),...
    'PlineLocation',round(abs(maxlat-minlat)/5),...%Draw a latitude line at the specified value
    'MLabelParallel','south' ...
    );

%Draw land
land=shaperead('landareas.shp','UseGeoCoords',true);
geoshow(ax,land,'FaceColor',[0.5 0.7 0.5])
h=geoshow(LL{1,sub}, R, 'DisplayType', 'texturemap');%Display geographic data
h.AlphaDataMapping = 'none'; % interpet alpha values as transparency values
h.FaceAlpha = 'texturemap'; % Indicate that the transparency can be different each pixel
alpha(h,double(~isnan(LL{1,sub})))
set(gca,'box','on')
title(['DIXSG Hour','-',num2str(sub)]);
colormap(ax1,parula(lev+1)); % Generate colors based on the value of lev
cb = colorbar;
clim([0, lev+1]); %Set the color bar range
cb.Location = 'eastoutside'; %Set the position of the color bar
cb.TickLength = 0.02; %Set the length of the color bar scale
cb.Ticks = linspace(0.5, lev+0.5, lev+1);
cb.TickLabels = cellstr(num2str((0:lev)', '%d')); %Set color bar scale labels
end