function Grid = Get_grid(IPP_pair, FDIXSG, minlon, minlat, maxlon, maxlat, numlon, numlat, dlon, dlat, H, g)
%% Divide the grid according to intervals and step sizes
%Divide the data into regular grids according to the center latitude and longitude
% INPUT:
%     IPP_pair: Center point latitude and longitude
%     FDIXSG: sensitivity
%     minlon, minlat, maxlon, maxlat,dlon, dlat,:Parameters for dividing the grid
%     numlon, numlat:Number of divisions based on longitude and latitude
%     H:field name
%     g:Number of data columns
% OUTPUT:
%     Grid:Divided data grid
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
Grid = cell(1, 24);
lon_edges = minlon:dlon:maxlon;
lat_edges = minlat:dlat:maxlat;
countnan=0;
for k = 1:24 % 24 Hours
    LL = cell(1, numlon); %Store data within each longitude range
    for lon_index = 1:numlon %Longitude index
        data = struct(); %Initialize the structure array to store data within each latitude range
        for lat_index = 1:numlat %Latitude index
            data.(H{1, lat_index}) = []; %Initialize each field in the structure array as an empty array
        end

        for i = 1:g %Column
            for j = 120*k-(120-1):120*k %Hourly data
                if IPP_pair.lon(j,i)~=0 && IPP_pair.lat(j,i)~=0 && ~isnan(FDIXSG(j,i)) && ~isnan(IPP_pair.lon(j,i)) && ~isnan(IPP_pair.lat(j,i)) %The current data is within the longitude range
                    if IPP_pair.lon(j,i) > lon_edges(lon_index) && IPP_pair.lon(j,i) <= lon_edges(lon_index + 1)
                        for lat_index=1:numlat
                            if IPP_pair.lat(j,i) > lat_edges(lat_index) && IPP_pair.lat(j,i) <= lat_edges(lat_index+1)
                                data.(H{1, lat_index})(end+1) = FDIXSG(j,i); %Store data in corresponding fields
                            end
                        end
                    end
                end
            end
        end
        %Delete null value fields in the structure
        empty_fields = cellfun(@isempty, struct2cell(data));
        data = rmfield(data, H(empty_fields));
        if ~isempty(fieldnames(data)) %Check if the structure array is empty
            LL{1, lon_index} = data; %Store data in LL
        end
    end
    Grid{1, k} = LL; %All latitude and longitude intervals within 1 hour for a total of 24 hours
    if all(cellfun(@isempty, LL))
        countnan=countnan+1;
    end
end
if countnan==24
    Grid=[];
end
end
