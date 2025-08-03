function [data] = preprocess(data,options)
%% Main GNSS Observation Data Preprocessing Function
% Function:
%     This function serves as the core pipeline for preprocessing GNSS
%     observation data. It orchestrates a sequence of sub-functions in a
%     standard and logical order to clean, compute, filter, and refine the
%     input data. The entire workflow is designed to prepare high-quality
%     observations for subsequent precise positioning algorithms (e.g., PPP).
%
%     The processing steps include:
%     1.  `elm_badclk`:   Eliminates observations with bad satellite clock or
%                         coordinate data, performing an initial data cleaning.
%     2.  `decimation`:   Decimates the data according to the desired interval
%                         or removes data outside the specified processing time window.
%     3.  `cal_sat`:      Calculates the precise satellite positions for each epoch.
%     4.  `elv_mask`:     Applies an elevation mask to remove observations from
%                         satellites below a defined cutoff angle.
%     5.  `cs_detect`:    Performs cycle slip detection and correction, a critical
%                         step for carrier phase processing (e.g., using GF or
%                         MW combinations).
%     6.  `clk_jmp2`:     (Optional) Detects and handles receiver clock jumps.
%     7.  `outlier`:      Detects and removes outliers (gross errors) from the
%                         observation data.
%     8.  `smoothing`:    (Optional) Performs carrier-smoothing of pseudorange
%                         measurements, using the precise carrier phase to reduce
%                         noise in the code measurements (also known as the Hatch filter).
%
% INPUT:
%     data:       A structure containing the raw or semi-processed GNSS
%                 observation data.
%     options:    A configuration structure that controls the behavior of the
%                 preprocessing steps. For example:
%                 - options.elevation_mask: The elevation cutoff angle.
%                 - options.cs_method:      The chosen cycle slip detection method.
%                 - options.clkjump:        Flag to enable/disable receiver clock
%                                           jump detection (1 for enable, 0 for disable).
%                 - options.codsmth:        Flag to enable/disable code smoothing
%                                           (1 for enable, 0 for disable).
%
% OUTPUT:
%     data:       The fully preprocessed data structure, ready for subsequent
%                 positioning calculations.
%
%% ---------------------------------------------------------------------

[data]  = elm_badclk(data);%Clear observation data with satellite coordinate errors

[data] = decimation(data,options);%Clear observation data that is not within the calculation range

[data] = cal_sat(data);%Calculate the actual satellite coordinates

[data] = elv_mask(data,options);%Cut off angle

[data] = cs_detect(data,options);%Cycle jump detection and repair, optional GF method or MW method

if options.clkjump == 1
    [data] = clk_jmp2(data);%Receiver clock jump
end

[data] = outlier(data);%Exclude observation data with gross errors

if options.codsmth == 1
    [data] = smoothing(data);%Carrier pseudorange smoothing
end
end

