function [data] = decimation(data,options)
%% Clip Observation Data to a Specific Time Window
% Function:
%     This function truncates the GNSS observation data based on a user-defined
%     start time (`from`) and end time (`to`) specified in the options structure.
%     It removes all observation epochs that fall outside this designated time
%     window, allowing the subsequent analysis and processing to focus only on
%     the period of interest.
%
% INPUT:
%     data:       A structure containing GNSS observation data, with matrices
%                 like ep, p1, p2, l1, l2, st under the `data.obs` field.
%     options:    A configuration structure that must contain the following fields:
%                 - options.from: The desired start time for processing (Second of Day).
%                 - options.to:   The desired end time for processing (Second of Day).
%
% OUTPUT:
%     data:       The modified data structure, where the observation data now
%                 contains only the epochs within the [from, to] time interval.
%

if (options.from<data.obs.ep(1,1)) || (options.to>data.obs.ep(end,1))
    errordlg('Observation file may not contain data between the chosen interval.','Epoch Interval Error');
    error   ('Observation file may not contain data between the chosen interval.');
end

f = (options.from - data.obs.ep(1,1));
l = (options.to - data.obs.ep(1,1));
fe = round(f/data.inf.time.int) + 1; %first epoch
le = round(l/data.inf.time.int) + 1; %last epoch

if le<size(data.obs.st,1)
    data.obs.p1((le+1):end,:) = [];
    data.obs.p2((le+1):end,:) = [];
    data.obs.l1((le+1):end,:) = [];
    data.obs.l2((le+1):end,:) = [];
    data.obs.ep((le+1):end,:) = [];
    data.obs.st((le+1):end,:) = [];
end

if fe~=1
    data.obs.p1(1:(fe-1),:) = [];
    data.obs.p2(1:(fe-1),:) = [];
    data.obs.l1(1:(fe-1),:) = [];
    data.obs.l2(1:(fe-1),:) = [];
    data.obs.ep(1:(fe-1),:) = [];
    data.obs.st(1:(fe-1),:) = [];
end
end

