function data = checkAnomalies(data)
%% checkAnomalies - Detects anomalies and removes invalid data arcs.
% INPUT:
%     data: A [numRows x numCols] matrix where each column is a time series
%           which may contain anomalies.
% OUTPUT:
%     data: The cleaned data matrix with anomalies and short arcs replaced by NaN.
%
% written by Zhang P. et al., 2024/08
%% ----------------------- ------------------------------------

DIFF_THRESHOLD = 2.5;

MIN_ARC_LENGTH = 10; 

if isempty(data)
    return;
end

[numRows, numCols] = size(data);

A_with_nan = [nan(1, numCols); data];
row_differences = diff(A_with_nan);
final_result = abs(row_differences);

anomaly_indices = find(final_result > DIFF_THRESHOLD);

if ~isempty(anomaly_indices)

    for i = 1:length(anomaly_indices)
        [row, col] = ind2sub(size(final_result), anomaly_indices(i));
        
        if isnan(data(row, col))
            continue; 
        end

        arc_end_row = row;
        while arc_end_row + 1 <= numRows && ~isnan(data(arc_end_row + 1, col))
            arc_end_row = arc_end_row + 1;
        end

        data(row:arc_end_row, col) = NaN;
    end
end


%% ========================================================================
%  STEP 2: SHORT ARC REMOVAL
% =========================================================================

% Now, iterate through the data column by column to find and remove short arcs.
for col = 1:numCols
    
    row_idx = 1;
    while row_idx <= numRows
        
        % --- Skip over existing NaNs to find the start of a potential arc ---
        if isnan(data(row_idx, col))
            row_idx = row_idx + 1;
            continue;
        end
        
        % If we are here, we found the start of an arc
        arc_start_row = row_idx;
        
        % --- Find the end of this arc ---
        arc_end_row = arc_start_row;
        while arc_end_row + 1 <= numRows && ~isnan(data(arc_end_row + 1, col))
            arc_end_row = arc_end_row + 1;
        end
        
        % --- Check the arc's length ---
        arc_length = arc_end_row - arc_start_row + 1;
        
        % If the arc is shorter than the minimum required length, remove it
        if arc_length < MIN_ARC_LENGTH
            data(arc_start_row:arc_end_row, col) = NaN;
        end
        
        % --- Move the main pointer to the position after this arc ---
        % This is crucial for efficiency and correctness.
        row_idx = arc_end_row + 1;
        
    end % end while loop for rows
end % end for loop for columns

end % end function
