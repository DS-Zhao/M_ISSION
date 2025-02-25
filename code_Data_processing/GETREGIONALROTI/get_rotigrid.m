function grid = get_rotigrid(lats, lons, range, ROTI, grid,num_lats,num_lons)
%% Update the ROTI value of the grid: take the average
% INPUT:
%     lats,lons: The latitude and longitude of the corresponding puncture point for that hour
%     range: Grid division setting
%            Includes: minlat, maxlat, minlon, maxlon, dlat, dlon
%     ROTI：ROTI data for that hour
%     grid：Initialize the grid
%     num_lats,num_lons,: Number of latitude and longitude divisions
% OUTPUT:
%     grid：Grid after data update
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
lat_min = range(1);
lon_min = range(3);
lat_resolution = range(5);
lon_resolution = range(6);
for i = 1:size(ROTI, 2)
    for j = 1:size(ROTI, 1)
        if ~isnan(ROTI(j, i))
            lat_index = ceil((lats(j, i) - lat_min) / lat_resolution) ;
            lon_index = ceil((lons(j, i) - lon_min) / lon_resolution) ;

            if 0<lat_index&&lat_index<= num_lats && 0<lon_index&&lon_index <= num_lons
                if isnan(grid(lat_index, lon_index))
                    grid(lat_index, lon_index) = ROTI(j, i);
                else
                    grid(lat_index, lon_index) = (grid(lat_index, lon_index) + ROTI(j, i)) / 2;
                end
            end
        end
    end
end
end