function  corrected_data=model_cor(data,model,folderPath)
%% Observing data for model correction
% INPUT:
%     data: Observation data obtained by calling PPPH
%     model: Call PPPH to obtain the model correction value
%     folderPath：folder path
% OUTPUT:
%     corrected_data:Observation data after model correction
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
[error1, error2, error3, error4,len1, len2, len3, len4] = extract_errors(model);
corrected_data = correct_observations(data, error1, error2, error3, error4);
corrected_data = correct_observations(corrected_data, len1, len2, len3, len4);
grandPath=fileparts(fileparts(folderPath));
year=num2str(data.inf.time.first(1));
doy=data.inf.time.doy;
[~, fileName, ~] = fileparts(data.files.rinex);
name=['COR_',fileName(1:4),year(3:4),num2str(doy,'%03d'),'.mat'];
newDirectoryName = ['raw_OBS/',year(3:4),num2str(doy,'%03d')];
cordata_path= fullfile(grandPath, newDirectoryName);
if exist(cordata_path,'dir')==0 
    mkdir(cordata_path);
end  
save([cordata_path,'\',name],'corrected_data','-mat');
end
function [error1, error2, error3, error4,len1, len2, len3, len4] = extract_errors(model)
    %Retrieve the number of rows for model data
    numRows = size(model, 1);
    
    %Initialize four error matrices
    error1 = [];
    error2 = [];
    error3 = [];
    error4 = [];
    len1=[]; len2=[]; len3=[]; len4=[];
    
    %Initialize the counter for the number of occurrences and distinguish which band is being corrected
    machineCount = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    
    % Read model data line by line
    for i = 1:numRows
        machineNumber = model(i, 4);  
        rowNumber = floor(model(i, 3) / 30) + 1;  
        errorValue = model(i, 7);  
        lenValue=model(i, 14);
       
        if isKey(machineCount, machineNumber)
            count = machineCount(machineNumber);
            if count == 4
                machineCount(machineNumber) = 1;
            else
                machineCount(machineNumber) = count + 1;
            end
        else
            machineCount(machineNumber) = 1;
        end
        
        %Determine which variable to store based on the number of occurrences
        count = machineCount(machineNumber);
        switch count
            case 1
                error1(rowNumber, machineNumber) = errorValue;
                len1(rowNumber, machineNumber)=lenValue;
            case 2
                error2(rowNumber, machineNumber) = errorValue;
                len2(rowNumber, machineNumber)=lenValue;
            case 3
                error3(rowNumber, machineNumber) = errorValue;
                len3(rowNumber, machineNumber)=lenValue;
            case 4
                error4(rowNumber, machineNumber) = errorValue;
                len4(rowNumber, machineNumber)=lenValue;
            otherwise
                error('Unexpected error: machine number %d has an invalid count %d', machineNumber, count);
        end
    end
end

function corrected_data = correct_observations(data, error1, error2, error3, error4)
    %Obtain the size of observation data
    [numObsRows, numObsCols] = size(data.obs.p1);
    
    %Obtain the size of error data
    [numErrRows1, numErrCols1] = size(error1);
    [numErrRows2, numErrCols2] = size(error2);
    [numErrRows3, numErrCols3] = size(error3);
    [numErrRows4, numErrCols4] = size(error4);
    
    %Initialize the corrected data
    data.obs.p1(data.obs.p1==0)=nan;
    data.obs.p2(data.obs.p2==0)=nan;
    data.obs.l1(data.obs.l1==0)=nan;
    data.obs.l2(data.obs.l2==0)=nan;
    corrected_data = data;
    %Traverse each observation matrix and corresponding error matrix
    for row = 1:numObsRows
        for col = 1:numObsCols
            %Check Err1 and make corrections
            if row > numErrRows1 || col > numErrCols1 || error1(row, col) == 0 || isnan(error1(row, col))
                corrected_data.obs.p1(row, col) = NaN;
            else
                corrected_data.obs.p1(row, col) = data.obs.p1(row, col) - error1(row, col);
            end
            
            %Check error 2 and make corrections
            if row > numErrRows2 || col > numErrCols2 || error2(row, col) == 0 || isnan(error2(row, col))
                corrected_data.obs.p2(row, col) = NaN;
            else
                corrected_data.obs.p2(row, col) = data.obs.p2(row, col) - error2(row, col);
            end
            
            %Check error 3 and make corrections
            if row > numErrRows3 || col > numErrCols3 || error3(row, col) == 0 || isnan(error3(row, col))
                corrected_data.obs.l1(row, col) = NaN;
            else
                corrected_data.obs.l1(row, col) = data.obs.l1(row, col) - error3(row, col);
            end
            
            %Check Err4 and make corrections
            if row > numErrRows4 || col > numErrCols4 || error4(row, col) == 0 || isnan(error4(row, col))
                corrected_data.obs.l2(row, col) = NaN;
            else
                corrected_data.obs.l2(row, col) = data.obs.l2(row, col) - error4(row, col);
            end
        end
    end
end
