function [atx] = r_antx(f_anten,inf,options)
%% Read ANTEX Antenna Phase Center File
% Function:
%     This function parses an ANTEX format file to extract antenna Phase
%     Center Offset (PCO) values for both satellites and a specific receiver.
%
%
% INPUT:
%     f_anten:  The full path to the ANTEX file.
%     inf:      An information structure. It must contain the field `inf.ant.type`,
%               which is a string specifying the name of the receiver antenna
%               type (usually including radome) to look up in the file.
%     options:  An options structure. The `options.system` substructure is used
%               to specify which satellite systems to process (e.g.,
%               `options.system.gps = 1`).
%
% OUTPUT:
%     atx:      A structure containing the antenna PCO data.
%               - atx.sat.neu: Satellite antenna PCOs, stored in a 3D matrix
%               - atx.rcv.neu: Receiver antenna PCOs, stored in a 3D matrix 
%
%% ---------------------------------------------------------------------

[fid,errmsg] = fopen(f_anten);

if any(errmsg)
    errordlg('Antex file can not be opened.','Antex file error');
    error   ('Antex file error');
end

type = inf.ant.type;
sn   = 105;
sat_neu = NaN(sn,3,2);
rcv_neu = NaN( 1,3,4);

linenum = 0;
while ~feof(fid)
    
    tline = fgetl(fid);
    linenum = linenum + 1;
    tag   = strtrim(tline(61:end));
    
    if strcmp(tag,'START OF ANTENNA')
        tline = fgetl(fid);
        linenum = linenum + 1;
        tag   = strtrim(tline(61:end));
        
        if strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'G') && (options.system.gps == 1)
            sat_no = sscanf(tline(22:23),'%d');
            if sat_no>32
                continue
            end
            while ~strcmp(tag,'END OF ANTENNA')
                tline = fgetl(fid);
                linenum = linenum + 1;
                tag   = strtrim(tline(61:end));
                if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G01')
                    frq_no = 1; %L1
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G02')
                    frq_no = 2; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                end
            end
            
        elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'R') && (options.system.glo == 1)
            sat_no = 32 + sscanf(tline(22:23),'%d');
            if sat_no>58 
                continue
            end
            while ~strcmp(tag,'END OF ANTENNA')
                tline = fgetl(fid);
                linenum = linenum + 1;
                tag   = strtrim(tline(61:end));
                if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R01')
                    frq_no = 1; %L1
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R02')
                    frq_no = 2; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                end
            end
            
        elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'E') && (options.system.gal == 1)
            sat_no = 58 + sscanf(tline(22:23),'%d');
            if sat_no>88
                continue
            end
            while ~strcmp(tag,'END OF ANTENNA')
                tline = fgetl(fid);
                linenum = linenum + 1;
                tag   = strtrim(tline(61:end));
                if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'E01')
                    frq_no = 1; %L1
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'E05')
                    frq_no = 2; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                end
            end
            
        elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'C') && (options.system.bds == 1)
            sat_no = 88 + sscanf(tline(22:23),'%d');
            if sat_no>105
                continue
            end
            while ~strcmp(tag,'END OF ANTENNA')
                tline = fgetl(fid);
                linenum = linenum + 1;
                tag   = strtrim(tline(61:end));
                if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'C01')
                    frq_no = 1; %L1
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'C07')
                    frq_no = 2; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                end
            end
            
        elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(strtrim(tline(1:20)),type)
            while ~strcmp(tag,'END OF ANTENNA')
                tline = fgetl(fid);
                linenum = linenum + 1;
                tag   = strtrim(tline(61:end));
                if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G01')
                    frq_no = 1; %L1
                    tline = fgetl(fid);
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        rcv_neu(1,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G02')
                    frq_no = 2; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        rcv_neu(1,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R01')
                    frq_no = 3; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        rcv_neu(1,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R02')
                    frq_no = 4; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        rcv_neu(1,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                end
            end 
        end
    else
        continue
    end
end

atx.sat.neu = sat_neu./1000;
atx.rcv.neu = rcv_neu./1000;

fclose('all');
end