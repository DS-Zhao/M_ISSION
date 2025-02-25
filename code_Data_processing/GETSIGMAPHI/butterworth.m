
function filtered_data=butterworth(corrected_data,cutoff_freq)
%% 6th butterworth filter
% INPUT:
%     corrected_data: Observation data after model correction
%     cutoff_freq: Filter threshold
% OUTPUT:
%     filtered_data:Filtered data
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
order = 6; %Filter order
%Design Butterworth filter
[b, a] = butter(order, cutoff_freq,"high");


%Initialize the filtered data structure
filtered_data = corrected_data;

%Traverse each observation matrix and filter the data
fields = {'l1', 'l2'};
[~,wavl] = frequencies;
for i = 1:numel(fields)
    field = fields{i};
    data2 = corrected_data.obs.(field);
    data2 = data2(:,1:105);
    data2=2*pi*data2./wavl(:,i)';
    %Check if the data is a column vector, otherwise transpose it to a column vector
    if isrow(data2)
        data2 = data2(:,1:105)';
    end

    %Initialize filtered data
    filtered_col = NaN(size(data2));

    %Find continuous non empty and non infinite data segments
    valid_idx = isfinite(data2);
    qq=zeros(1,size(valid_idx,2));
    d = diff([qq;valid_idx; qq]);
    start_idx = find(d == 1);
    idx = find(d==-1)-1;
    end_idx = find(d==-1)-floor(idx / 2881)-1;


    %Filter each continuous non empty and non infinite data segment
    for j = 1:length(start_idx)
        segment = data2(start_idx(j):end_idx(j));
        if length(segment) > 18 %The length of the data segment needs to be greater than the order of the filter
            filtered_segment = filtfilt(b, a, segment);
            filtered_col(start_idx(j):end_idx(j)) = filtered_segment;
        else
            filtered_col(start_idx(j):end_idx(j)) = nan; %Keep it as it is
        end
    end
    %Update filtered data
    filtered_data.obs.(field) = filtered_col;

end

