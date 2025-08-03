function [ obs,inf ] = read_obsf(files,options)
%% Read RINEX Observation File
% Function:
%     This function reads a RINEX observation (OBS) file. It automatically 
%     detects the RINEX version (supports 2.x and 3.x) and calls the 
%     appropriate sub-functions to parse the header and observation data. 
%
% INPUT:
%     files:    A structure containing the input file paths.
%               - files.rinex: Full path to the RINEX observation file.
%     options:  A structure containing processing options.
%               - options.dcb: Flag to enable/disable DCB correction (1: enable, 0: disable).
%
% OUTPUT:
%     obs:      A structure containing the observation data read from the RINEX file.
%     inf:      A structure containing various information parsed from the RINEX file header.
%% ---------------------------------------------------------------------

[ver] = r_rnxvers(files.rinex);

if ver>=3
    
    [inf] = r_rnxheadv3(files.rinex);
    
    
    [obs] = r_rnxobsv3(files.rinex,inf,options);

elseif ver>=2
    
    [inf] = r_rnxheadv2(files.rinex);
    
    
    [obs] = r_rnxobsv2(files.rinex,inf,options)   ;
   
else
    errordlg('RINEX version is not valid !','RINEX version error');
    error('RINEX version is not valid !');
end 
end

