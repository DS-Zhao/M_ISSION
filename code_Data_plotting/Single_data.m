function SS = Single_data(full_path)
%% Filter file types to match drawing functions
% This function processes a data file and extracts relevant information based on
% INPUT:
%     full_path: The complete path of the data file to be processed
% OUTPUT：
%     SS：Effective satellite
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
    A = load(full_path);
    A = struct2cell(A);
    A = A{1};
    
    [~, fileName, ~] = fileparts(full_path);
    
    last_four = fileName(end-3:end);% Get the last four characters of the file name
    
    % Determine the type of data based on the file suffix
    if strcmp('IXSG', last_four)
        SS = "DIXSG";
    else
        % Handle 'aphi' suffix files
        if strcmp('aphi', last_four)
            % Check if the struct contains 'L1' or 'L2' field and extract accordingly
            if ismember('L1', fieldnames(A))
                A = A.L1;
            else
                A = A.L2;
            end
        end
        
        % Call shiyan1 function to find columns that are not entirely NaN
        empty_cols = shiyan1(A);
        % Call shiyan2 function to simplify the data
        simplified = shiyan2(empty_cols);
        
        % Initialize SS as an empty string
        SS = '';
        % Iterate through the simplified data and concatenate it into SS
        for i = 1:numel(simplified)
            % Convert cell content to string
            if ischar(simplified{i})
                SS = [SS, simplified{i}, ' '];
            elseif isnumeric(simplified{i})
                % Handle numeric data by enclosing in brackets and separating with commas
                if isvector(simplified{i}) && length(simplified{i}) == 2
                    SS = [SS, '[', num2str(simplified{i}(1)), ',', num2str(simplified{i}(2)), ']', ','];
                else
                    SS = [SS,'[', num2str(simplified{i}), ']', ','];
                end
            end
        end
        % Remove the trailing comma and space
        if ~isempty(SS)
            SS = SS(1:end-1);
        end
    end
end


function empty_cols = shiyan1(A)
%% Helper function to find columns that are not entirely NaN
    empty_cols = [];
    for col = 1:size(A, 2)
        if ~all(isnan(A(:, col)))
            empty_cols = [empty_cols, col];
        end
    end
end


function simplified = shiyan2(arr)
%% Helper function to simplify the data
    simplified = {};
    start_idx = 1;
    end_idx = 1;
    for i = 2:length(arr)
        if arr(i) - arr(i-1) ~= 1
            % Check if the range is valid and add it to the simplified array
            if end_idx - start_idx >= 2
                simplified{end+1} = [arr(start_idx), arr(end_idx)];
            else
                simplified{end+1} = arr(start_idx:end_idx);
            end
            % Update the start and end indices for the next range
            start_idx = i;
            end_idx = i;
        else
            end_idx = i;
        end
    end
    % Add the last range if it's valid
    if end_idx - start_idx >= 2
        simplified{end+1} = [arr(start_idx), arr(end_idx)];
    else
        simplified{end+1} = arr(start_idx:end_idx);
    end
end