function data = checkAnomalies(data)
%% Input the data and replace anomalies with NaN
% INPUT:
%     data: input data which may contain anomalies
% OUTPUT:
%     data: data with anomalies replaced by NaN

%% written by Zhang P. et al., 2024/08
%% ---------------------------------------------------------------------

% Get the size of the data
[numRows, numCols] = size(data);

% Initialize thresholds and window size
threshold_increase = 2;
threshold_decrease = 2;
window_size = 10;

% Check data column by column
for col = 1:numCols
    validIndices = ~isnan(data(:, col));

    if all(~validIndices)
        continue;
    end

    processed = false(numRows, 1);

    for i = 1:numRows
        if processed(i) || isnan(data(i, col))
            continue;
        end

        % Find the start of a new segment
        if i == 1||isnan(data(i-1, col))
            startIdx = i;
        else
            startIdx = i - 1;
        end

        % Find the end of the segment
        endIdx = i;
        while endIdx < numRows && ~isnan(data(endIdx+1, col)) && ...
               (data(endIdx+1, col) - data(endIdx, col) <= threshold_increase)
            endIdx = endIdx + 1;
        end

        % Check if the segment is within the window size
        if endIdx - startIdx + 1 <= window_size
            continue;
        end

        % Check for a decrease within the window size after the increase
        foundDecrease = false;
        for j = endIdx:-1:startIdx + 1
            if data(j, col) - data(j-1, col) < -threshold_decrease
                foundDecrease = true;
                break;
            end
            if j == startIdx + 1
                break;
            end
        end

        % If a decrease is found, replace the segment with NaN
        if foundDecrease
            data(startIdx:endIdx, col) = NaN;
            processed(startIdx:endIdx) = true;
        end
    end
end
end