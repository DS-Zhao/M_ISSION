function [data] = elm_badclk(data)
%% Eliminate Observations with Bad Satellite Clock Data
% Function:
%     This function identifies and flags GNSS observations that correspond
%     to invalid or missing satellite clock corrections. In precise positioning,
%     accurate satellite clock data is essential. Any epoch lacking a valid
%     clock correction renders the satellite's observation at that epoch unusable.
%
% INPUT:
%     data:       A GNSS data structure containing the following fields:
%                 - data.clk: [num_epochs x num_sats] matrix of satellite clock corrections.
%                 - data.obs.st: [num_obs_epochs x num_sats] status matrix; 1 for valid, 0 for invalid.
%                 - data.obs.ep: [num_obs_epochs x ...] matrix, where the first column is the time of observation (Second of Day).
%                 - data.inf.time.clkint: The time interval (in seconds) of the satellite clock product.
%
% OUTPUT:
%     data:       The updated data structure where the `data.obs.st` status matrix
%                 has been modified. Observations corresponding to bad clock data
%                 are now flagged with 0.
%
%% ---------------------------------------------------------------------
intr = data.inf.time.clkint;
for i=1:size(data.clk,2)
    if (any(~isnan(data.clk(:,i)))) && (any(isnan(data.clk(:,i)) | data.clk(:,i)==999999.999999))
        loc = find(isnan(data.clk(:,i)));
        for t=loc'
            sod = (t-1)*intr;
            ep  = data.obs.ep(:,1)==sod;
            data.obs.st(ep,i) = 0;
        end
    elseif any(isnan(data.clk(:,i)))
        data.obs.st(:,i) = 0;
    end
end

end

