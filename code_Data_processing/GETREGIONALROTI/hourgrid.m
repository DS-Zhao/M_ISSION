function ROTIgrid=hourgrid(L,B,ROTI,num_lats,num_lons,range,ROTIgrid)
%% Regional ROTI grid division
% INPUT:
%     L,B: Longitude and latitude of puncture point
%     ROTI：ROTI data
%     num_lats,num_lons,: Number of latitude and longitude divisions
%     range: Grid division setting
%            Includes: minlat, maxlat, minlon, maxlon, dlat, dlon
% OUTPUT:
%     ROTIgrid：Grid filled with ROTI data
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
for k = 1:24
    startIdx = 120 * (k - 1) + 1;
    endIdx = 120 * k;
    %Initialize an empty grid for the current ROTI file
    grid = nan(num_lats, num_lons);

    ROTI_segment = ROTI(startIdx:endIdx, :);
    lats_segment = B(startIdx:endIdx, :);
    lons_segment = L(startIdx:endIdx, :);

    %Update Grid
    grid=get_rotigrid(lats_segment, lons_segment, range, ROTI_segment, grid,num_lats,num_lons);


    if ~isempty(ROTIgrid{k})
        %Find the positions of non NaN and NAN
        validIndices = ~isnan(ROTIgrid{k}) & ~isnan(grid);
        nanIndices = isnan(ROTIgrid{k});
        %Calculate the average of non NaN values and update ROTigrid {k}
         ROTIgrid{k}(validIndices) = (ROTIgrid{k}(validIndices) + grid(validIndices)) / 2;
        % ROTIgrid{k}(validIndices) = max(ROTIgrid{k}(validIndices),grid(validIndices));
        ROTIgrid{k}(nanIndices & ~isnan(grid)) = grid(nanIndices & ~isnan(grid));
    else
        ROTIgrid{k} = grid;
    end
end
end