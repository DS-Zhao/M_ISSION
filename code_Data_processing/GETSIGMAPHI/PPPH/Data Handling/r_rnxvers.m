function [ver] = r_rnxvers(f_obs)
%% Read RINEX Version Number
% Function:
%     Opens the specified RINEX observation file and reads the
%     version number from its header.
%
% INPUT:
%     f_obs:    Path to the RINEX observation file.
%
% OUTPUT:
%     ver:      The version number of the RINEX file (e.g., 2.11, 3.03).
%% ---------------------------------------------------------------------

[fid,errmsg] = fopen(f_obs);

if any(errmsg)
    errordlg('OBSERVATION file can not be opened !','Observation File Error');
    error   ('OBSERVATION file can not be opened !');
end

while 1
    
    tline = fgetl(fid);
    tag  = strtrim(tline(61:end));
    if strcmp(tag,'RINEX VERSION / TYPE')
        ver = sscanf(tline(1:20),'%f');
        break
    end
end
fclose('all');
end