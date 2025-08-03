function [obs] = r_rnxobsv3(f_obs,inf,options,dcb)
%% Read RINEX Version 3.x Observation Data Body
% Function:
%     This function reads the data body of a RINEX version 3.x observation
%     (OBS) file. It processes the file epoch by epoch, using the header
%     information from `r_rnxheadv3` (especially the observation sequence
%     in `inf.seq`) to correctly extract dual-frequency pseudorange and
%     carrier phase observations for each satellite system.
%
%
% INPUT:
%     f_obs:    Path to the RINEX v3 observation file.
%     inf:      Header information structure, returned by `r_rnxheadv3`.
%     options:  A structure with processing options.
%               - options.system: Specifies systems to process (e.g., .gps=1, .glo=1).
%
% OUTPUT:
%     obs:      A structure containing all extracted and processed observation data.
%               All data matrices are of size [number of epochs x 105].
%               - obs.p1: P1/C1 pseudorange matrix (meters).
%               - obs.p2: P2/C2 pseudorange matrix (meters).
%               - obs.l1: L1 carrier phase matrix (meters).
%               - obs.l2: L2 carrier phase matrix (meters).
%               - obs.ep: Vector of epoch times (seconds of day).
%               - obs.st: Status matrix (1 if data exists for the sat/epoch, 0 otherwise).
%% ---------------------------------------------------------------------

narginchk(3,4)
[fid,errmsg] = fopen(f_obs);

if any(errmsg)
    errordlg('OBSERVATION file can not be opened !','Observation file error');
    error   ('OBSERVATION file error');
end

sno = 105;

if inf.time.last == 86400
    fi = inf.time.first(4,1)*3600 + inf.time.first(5,1)*60 + inf.time.first(6,1);
    la = inf.time.last;
    max = (la-fi)/inf.time.int + 1;
else
    fi = inf.time.first(4,1)*3600 + inf.time.first(5,1)*60 + inf.time.first(6,1);
    la = inf.time.last (4,1)*3600 + inf.time.last (5,1)*60 + inf.time.last (6,1);
    max = (la-fi)/inf.time.int + 1;
end

p1s  = NaN(max,sno);
p2s  = NaN(max,sno);
l1s  = NaN(max,sno);
l2s  = NaN(max,sno);
eps  = NaN(max, 1);
st   = zeros(max,sno);
[~,wavl] = frequencies;

epno = 0; % epoch no
lno  = 0; % line no
while ~feof(fid)

    tline = fgetl(fid);
    lno   = lno + 1;
    if ~isempty(tline) && strcmp(tline(1),'>')
        nep  = sscanf(tline(3:end),'%f');
        if size(nep,1)>4
            epno = epno + 1;
            eps(epno,1) = nep(4)*3600 + nep(5)*60 + nep(6);
            for i=1:nep(8)

                tline = fgetl(fid);
                lno   = lno + 1;
                % GPS
                if strcmp(tline(1),'G') && (options.system.gps == 1)
                    k   = sscanf(tline(2:3),'%d');
                    if k>32
                        continue
                    end
                    ono = inf.nob.gps;
                    lso = NaN(ono,1);
                    nu  = 0;
                    for u=4:16:size(tline,2)
                        nu = nu + 1;
                        tls = sscanf(tline(u:u+13),'%f');
                        if ~isempty(tls)
                            lso(nu,1) = tls;
                        end
                    end

                    if options.dcb == 1
                        p1s(epno,k) = lso(inf.seq.gps(1)) + dcb(k,1);
                    else
                        p1s(epno,k) = lso(inf.seq.gps(1));
                    end
                    p2s(epno,k) = lso(inf.seq.gps(2));
                    l1s(epno,k) = lso(inf.seq.gps(3))*wavl(k,1);
                    l2s(epno,k) = lso(inf.seq.gps(4))*wavl(k,2);
                    % GLONASS
                elseif strcmp(tline(1),'R') && (options.system.glo == 1)
                    k   = 32 + sscanf(tline(2:3),'%d'); % satellite prn number
                    if k>58
                        continue
                    end
                    ono = inf.nob.glo;               % observation number
                    lso = NaN(ono,1);
                    nu  = 0;
                    for u=4:16:size(tline,2)
                        nu = nu + 1;
                        tls = sscanf(tline(u:u+13),'%f');
                        if ~isempty(tls)
                            lso(nu,1) = tls;
                        end
                    end
                    % store the data
                    p1s(epno,k) = lso(inf.seq.glo(1));
                    p2s(epno,k) = lso(inf.seq.glo(2));
                    l1s(epno,k) = lso(inf.seq.glo(3))*wavl(k,1);
                    l2s(epno,k) = lso(inf.seq.glo(4))*wavl(k,2);
                    % GALILEO
                elseif strcmp(tline(1),'E') && (options.system.gal == 1)
                    k   = 58 + sscanf(tline(2:3),'%d'); % satellite prn number
                    if k>88
                        continue
                    end
                    ono = inf.nob.gal;               % observation number
                    lso = NaN(ono,1);
                    nu  = 0;
                    for u=4:16:size(tline,2)
                        nu = nu + 1;
                        tls = sscanf(tline(u:u+13),'%f');
                        if ~isempty(tls)
                            lso(nu,1) = tls;
                        end
                    end
                    % store the data
                    p1s(epno,k) = lso(inf.seq.gal(1));
                    p2s(epno,k) = lso(inf.seq.gal(2));
                    l1s(epno,k) = lso(inf.seq.gal(3))*wavl(k,1);
                    l2s(epno,k) = lso(inf.seq.gal(4))*wavl(k,2);
                    % BEIDOU
                elseif strcmp(tline(1),'C') && (options.system.bds == 1)
                    k   = 88 + sscanf(tline(2:3),'%d'); % satellite prn number
                    if k>105
                        continue
                    end
                    ono = inf.nob.bds;               % observation number
                    lso = NaN(ono,1);
                    nu  = 0;
                    for u=4:16:size(tline,2)
                        nu = nu + 1;
                        tls = sscanf(tline(u:u+13),'%f');
                        if ~isempty(tls)
                            lso(nu,1) = tls;
                        end
                    end
                    % store the data
                    p1s(epno,k) = lso(inf.seq.bds(1));
                    p2s(epno,k) = lso(inf.seq.bds(2));
                    l1s(epno,k) = lso(inf.seq.bds(3))*wavl(k,1);
                    l2s(epno,k) = lso(inf.seq.bds(4))*wavl(k,2);
                end
            end
        end

    end
end

all = p1s + p2s + l1s + l2s;
st(~isnan(all)) = 1;

if max>epno
    p1s(epno+1:max,:) = [];
    p2s(epno+1:max,:) = [];
    l1s(epno+1:max,:) = [];
    l2s(epno+1:max,:) = [];
    eps(epno+1:max,:) = [];
    st(epno+1:max,:)  = [];
end

obs.p1 = p1s;
obs.p2 = p2s;
obs.l1 = l1s;
obs.l2 = l2s;
obs.ep = eps;
obs.st = st;

fclose('all');
end

