function [clk] = r_clck(f_clock,options)
%% Read RINEX Precise Clock File
% Function:
%     This function parses a RINEX format precise clock (.clk) file.
%     Based on the specified system options, it reads and processes satellite
%     clock data for GPS, GLONASS, Galileo, and BDS.
%
% INPUT:
%     f_clock:  The full path to the precise clock (.clk) file.
%     options:  An options structure.
%               - options.system: A sub-structure specifying which satellite
%                                 systems to read (e.g., .gps=1, .glo=1).
%               - options.clck_int: The sampling interval of the clock file
%                                   in seconds, used to calculate the total
%                                   number of epochs and the row index for
%                                   each data point.
%
% OUTPUT:
%     clk:      A 2D matrix [epoch, satellite_index] storing the precise
%               clock bias values for all satellites, in seconds. Unfilled
%               elements will remain NaN.
%% ---------------------------------------------------------------------

[fid,errmsg] = fopen(f_clock);

if any(errmsg)
    errordlg('Clock file can not be opened.','Clock file error');
    error   ('Clock file error');
end

sn = 105;
tn = 86400/options.clck_int;%30s
clk = NaN(tn,sn);

while ~feof(fid)
    tline = fgetl(fid);
    
    if strcmp(tline(1:4),'AS G') && (options.system.gps == 1)
        new = sscanf(tline(5:end),'%f',[1,10]);%02  2017  7 10  0  0  0.000000  1    0.380593971661E-03
        sat_no = new(1);%2
        if sat_no<33
            epoch  = (new(5)*3600 + new(6)*60 + new(7))/options.clck_int + 1;
            clk(epoch,sat_no) = new(9);
        else
            continue
        end
    elseif strcmp(tline(1:4),'AS R') && (options.system.glo == 1)
        new = sscanf(tline(5:end),'%f',[1,10]);
        sat_no = 32 + new(1);
        if sat_no<59
            epoch  = (new(5)*3600 + new(6)*60 + new(7))/options.clck_int + 1;
            clk(epoch,sat_no) = new(9);
        else
            continue
        end
    elseif strcmp(tline(1:4),'AS E') && (options.system.gal == 1)
        new = sscanf(tline(5:end),'%f',[1,10]);
        sat_no = 58 + new(1);
        if sat_no<89
            epoch  = (new(5)*3600 + new(6)*60 + new(7))/options.clck_int + 1;
            clk(epoch,sat_no) = new(9);
        else
            continue
        end
    elseif strcmp(tline(1:4),'AS C') && (options.system.bds == 1)
        new = sscanf(tline(5:end),'%f',[1,10]);
        sat_no = 88 + new(1);
        if sat_no<106
            epoch  = (new(5)*3600 + new(6)*60 + new(7))/options.clck_int + 1;
            clk(epoch,sat_no) = new(9);
        else
            continue
        end
    end 
end
fclose('all');
end