function [validData,xValues]=time24h(data)
%% Convert horizontal axis time coordinates (1-24 hours) based on data length
%INPUT:
%     data: Data to be drawn
%OUTPUT:
%     validData：Valid data without null values
%     xValues：The x-axis time label corresponding to the valid data
%% written by Zhang P. et al., 2024/11
%% -----------------------------------------------------------------------
% Find the indexes of all non NaN rows
nonNaNRows = ~all(isnan(data), 2); 
rowIndices = find(nonNaNRows); %Retrieve the index of non NaN rows

%Check for non NaN rows
if ~isempty(rowIndices)
    %Retrieve the indexes of the first and last rows
    firstValidRow = rowIndices(1);
    lastValidRow = rowIndices(end);
    
end

timePerRow = 24 / size(data,1);  
   
validData = data(firstValidRow:lastValidRow, :);  
%Corresponding x time value (only within the valid range, taking into account the time interval)  
xValues = (firstValidRow:lastValidRow)' * timePerRow;%Transpose to column vector
end
  