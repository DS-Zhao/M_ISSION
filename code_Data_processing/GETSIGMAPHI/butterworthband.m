function filtered_data = butterworthband(corrected_data, band_freqs, fs)
%% 6th order Butterworth bandpass filter

% INPUT:
%   corrected_data: Structure containing observation data after model correction.
%                   It should have a field .obs with subfields .l1, .l2 etc.
%   band_freqs:     A 1x2 vector specifying the passband frequencies [f_low, f_high]
%                   in Hz. For example, [0.001, 0.015] for MSTIDs analysis with 30s data.
%   fs:             The sampling frequency of the data in Hz (e.g., 1/30 for 30s data).
%
% OUTPUT:
%   filtered_data:  Data structure with phase observations filtered.
%
%% written by Zhang P. et al., updated 2025/07
%% -----------------------------------------------------------------------
order = 6; % Filter order

nyquist_freq = fs / 2;
Wn = band_freqs / nyquist_freq; 
[b, a] = butter(order, Wn, 'bandpass'); 

filtered_data = corrected_data;

fields = {'l1', 'l2'}; 

for i = 1:numel(fields)
    field = fields{i};
    if ~isfield(corrected_data.obs, field)
        continue; % Skip if this field does not exist
    end

    data2 = corrected_data.obs.(field);
    filtered_cols = NaN(size(data2));

    % Process each column (satellite) independently
    for sat_idx = 1:size(data2, 2)
        sat_data = data2(:, sat_idx);

        % Find continuous non-NaN and non-infinite data segments
        valid_idx = isfinite(sat_data);
        d = diff([0; valid_idx; 0]); 
        start_indices = find(d == 1);
        end_indices = find(d == -1) - 1;

        % Filter each continuous segment
        for j = 1:length(start_indices)
            start_seg = start_indices(j);
            end_seg = end_indices(j);
            segment = sat_data(start_seg:end_seg);
            
            if length(segment) > 6 * order
                segment2=detrend(segment,3);%Third order polynomial detrended
                filtered_segment = filtfilt(b, a, segment2);
                filtered_cols(start_seg:end_seg, sat_idx) = filtered_segment;
            else
                filtered_cols(start_seg:end_seg, sat_idx)=nan;
            end
        end
    end
    % Update filtered data
    filtered_data.obs.(field) = filtered_cols;
end
end

