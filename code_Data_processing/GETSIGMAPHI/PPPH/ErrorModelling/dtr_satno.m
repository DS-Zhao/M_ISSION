function [satno] = dtr_satno(obs)
%% Determine the Number of Satellites per Epoch
% Function:
%     This function calculates the total number of satellites tracked at each 
%     observation epoch. 
% INPUT:
%     obs:        A structure containing GNSS observation data. This structure
%                 must include a field named 'st'. 
%                 'obs.st' is expected to be an m-by-n matrix, where:
%                 - m is the number of observation epochs.
%                 - n is the total number of satellites in the constellation (or max PRN).
%                 - The elements are typically 1 (if the satellite is tracked at 
%                   that epoch) and 0 (if it is not).
%
% OUTPUT:
%     satno:      An m-by-1 column vector, where satno(i) contains the total
%                 count of satellites observed at the i-th epoch.
%
%% ---------------------------------------------------------------------
m = size(obs.st,1);
satno = zeros(m,1);

for i=1:m
    satno(i,1) = sum(obs.st(i,:));
end

end

