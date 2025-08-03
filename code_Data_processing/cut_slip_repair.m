function [obs]=cut_slip_repair(path_obs,Sites_Info)
%% Remove small arc segments and cycle jump data
% INPUT:
%     path_obs: R_rinex:Observation data storage path
%     Sites_Info: name and coordinate information of the stations
% OUTPUT:
%     obs: Modified observation data
% SAVE:
%      */OBS_cut/doy/sitedoy.mat:Processed OBS data
%% modified by Zhang P. et al., 2024/08
%% *****************************************
doys=Sites_Info.doy;
stations=Sites_Info.name;
current_doy=num2str(unique(doys));
path_cutobs = fullfile(fileparts(path_obs), 'raw_OBS_cut');
list_obs=dir([path_cutobs,'/',current_doy,'/*.mat']);
list_obs=now_files(stations,list_obs);
len=length(list_obs);
coor=[];EA=[];
for i=1:len
    load([path_cutobs,'/',current_doy,'/',list_obs(i).name],'-mat');
    doy=list_obs(i).name(end-8:end-4);
    fields=fieldnames(obs);
    if ~isnan(find(strcmp(fields, 'GPSC1W' )))
        if ~(all(all(obs.GPSC1W==0)) || all(all(obs.GPSC2W==0)))
            obs=GPS_slip(obs,EA);
        end
    end
    if ~isnan(find(strcmp(fields, 'GLOC1P' )))
        if ~(all(all(obs.GLOC1P==0)) || all(all(obs.GLOC2P==0)))
            obs=GLO_slip(obs,EA);
        end
    end
    if ~isnan(find(strcmp(fields, 'BDSC2I')))
        if ~( all(all(obs.BDSC2I==0)) || all(all(obs.BDSC6I==0)))
            obs=BDS_slip(obs,EA);
        end
    elseif ~isnan(find(strcmp(fields, 'BDSC1I')))
        if ~( all(all(obs.BDSC1I==0)) || all(all(obs.BDSC7I==0)))
            obs=BDS1_slip(obs,EA);
        end
    end

    if ~isnan(find(strcmp(fields, 'GALC1C' )))
        if ~(all(all(obs.GALC1C==0)) || all(all(obs.GALC5Q==0)))
            obs=GAL_slip(obs,EA);
        end
    elseif ~isnan(find(strcmp(fields, 'GALC1X' )))
        if ~(all(all(obs.GALC1X==0)) || all(all(obs.GALC5X==0)))
            obs=GALX_slip(obs,EA);
        end
    end

    path_cutobs2 = fullfile(fileparts(path_obs), 'raw_OBS_cut',doy);
    filename1 = fullfile(path_cutobs2, list_obs(i).name);
    % Save data to file
    save(filename1, 'obs','coor','EA','-mat');
end
end
%%
function [obs] = GPS_slip(obs, EA)
%% GPS observations Cycle Slip Detection and Repair
% INPUT:
%      obs: struct of rinex observation files, containing fields like
%           GPSC1W, GPSC2W, GPSL1C, GPSL2W, and ideally a time vector.
%      EA:  struct containing elevation angles. Expected field: EA.GPSel
%           which should correspond to the GPS satellites.
%
% OUTPUT:
%      obs: Modified observation data with cycle slips repaired.

[epoch_num, sat_num] = size(obs.GPSC1W);
c = 299792458;                       
GPS_f1 = 1575.42e6;                   
GPS_f2 = 1227.60e6;                   
lamda1 = c / GPS_f1;                  
lamda2 = c / GPS_f2;                  
lamda_w = c / (GPS_f1 - GPS_f2);      
dt = 30; 
sig0 = sqrt(2 * (0.0027^2 + 0.0017^2)); 

iono_rate_per_hour = 4; 
dl = (iono_rate_per_hour) * (dt / 3600);

valid_mask = obs.GPSL1C ~= 0 & obs.GPSL2W ~= 0 & obs.GPSC1W ~= 0 & obs.GPSC2W ~= 0;

L_gf = zeros(epoch_num, sat_num);
L1_m = obs.GPSL1C .* lamda1;
L2_m = obs.GPSL2W .* lamda2;
L_gf(valid_mask) = L1_m(valid_mask) - L2_m(valid_mask);

N_wl = zeros(epoch_num, sat_num);
L_wl = (GPS_f1 * L1_m - GPS_f2 * L2_m) / (GPS_f1 - GPS_f2);
P_nl = (GPS_f1 * obs.GPSC1W + GPS_f2 * obs.GPSC2W) / (GPS_f1 + GPS_f2);
N_wl(valid_mask) = (L_wl(valid_mask) - P_nl(valid_mask)) / lamda_w;

for i = 1:sat_num 
    
    arc = Get_arc(valid_mask(:, i));
    if isempty(arc)
        continue;
    end
    
    arc( (arc(:, 2) - arc(:, 1) + 1) < 30, :) = [];
    if isempty(arc)
        continue;
    end
    
    for j = 1:size(arc, 1) 
        arc_start = arc(j, 1);
        arc_end   = arc(j, 2);

        arc_nwl_std = std(N_wl(arc_start:arc_end, i));
        
        for k = (arc_start + 1) : arc_end
            
            dmwc = 0; 
            dgfc = 0; 

            dgf = L_gf(k, i) - L_gf(k-1, i);
            elv = EA.GPSel(k, i);
            

            me = 1 + (10 * exp(-elv / 10));
            smg = sig0 * me;
            
            thresh_gf = (4 * smg) + dl;
            
            if abs(dgf) > thresh_gf
                dgfc = 1; 
            end

            window_start = max(arc_start, k - 30);

            if (k - window_start) >= 3
                past_nwl_data = N_wl(window_start : k-1, i);
                mmw = mean(past_nwl_data); 
                smw = std(past_nwl_data);  

                smw = max(smw, 0.01);
                
                dmw = N_wl(k, i) - mmw; 
                
                if abs(dmw) > (5 * smw)
                    dmwc = 1; 
                end
            end
            
            if dgfc || (dmwc && arc_nwl_std > 0.6)
                
                delta_Nwl = N_wl(k, i) - N_wl(k-1, i);
                delta_Lgf = L_gf(k, i) - L_gf(k-1, i);
                
                A = [1, -1; lamda1, -lamda2];
                L = [delta_Nwl; delta_Lgf];
                dN_float = A \ L;
                
                dN1_fix = round(dN_float(1));
                dN2_fix = round(dN_float(2));

                if dN1_fix ~= 0 || dN2_fix ~= 0
                    obs.GPSL1C(k:arc_end, i) = obs.GPSL1C(k:arc_end, i) - dN1_fix;
                    obs.GPSL2W(k:arc_end, i) = obs.GPSL2W(k:arc_end, i) - dN2_fix;
                    
                    L_gf(k:arc_end, i) = L_gf(k:arc_end, i) - (dN1_fix * lamda1 - dN2_fix * lamda2);
                    N_wl(k:arc_end, i) = N_wl(k:arc_end, i) - (dN1_fix - dN2_fix);
                end
                
            end
        end
    end
end
end


function [obs] = GLO_slip(obs, EA)
%% GLONASS observations Cycle Slip Detection and Repair
% INPUT:
%      obs: struct of rinex observation files, containing fields like
%           GLOC1P, GLOC2P, GLOL1P, GLOL2P.
%      EA:  struct containing elevation angles. Expected field: EA.GLOel.
%
% OUTPUT:
%      obs: Modified observation data with cycle slips repaired.

[epoch_num, sat_num] = size(obs.GLOC1P);
c = 299792458;

dt = 30; 
sig0 = sqrt(2 * (0.0027^2 + 0.0017^2)); 

iono_rate_per_hour = 4; 
dl = (iono_rate_per_hour) * (dt / 3600);

Fre = [1,-4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2, 0, 0];
if sat_num > length(Fre)
    warning('Number of GLONASS satellites exceeds the length of the frequency number array. Truncating satellite processing.');
    sat_num = length(Fre);
end

for i = 1:sat_num 

    k_freq = Fre(i); 
    f1 = (1602 + k_freq * 0.5625) * 1e6; 
    f2 = (1246 + k_freq * 0.4375) * 1e6; 
    lamda1 = c / f1; 
    lamda2 = c / f2; 
    lamda_w = c / (f1 - f2);   

    valid_mask = obs.GLOL1P(:, i) ~= 0 & obs.GLOL2P(:, i) ~= 0 & obs.GLOC1P(:, i) ~= 0 & obs.GLOC2P(:, i) ~= 0;

    L_gf = zeros(epoch_num, 1);
    N_wl = zeros(epoch_num, 1);

    if ~any(valid_mask)
        continue; 
    end

    L1_m = obs.GLOL1P(:, i) .* lamda1;
    L2_m = obs.GLOL2P(:, i) .* lamda2;
    L_gf(valid_mask) = L1_m(valid_mask) - L2_m(valid_mask);

    L_wl = (f1 * L1_m - f2 * L2_m) / (f1 - f2);
    P_nl = (f1 * obs.GLOC1P(:, i) + f2 * obs.GLOC2P(:, i)) / (f1 + f2);
    N_wl(valid_mask) = (L_wl(valid_mask) - P_nl(valid_mask)) / lamda_w;

    arc = Get_arc(valid_mask);
    if isempty(arc)
        continue;
    end

    arc( (arc(:, 2) - arc(:, 1) + 1) < 30, :) = [];
    if isempty(arc)
        continue;
    end

    for j = 1:size(arc, 1) 
        arc_start = arc(j, 1);
        arc_end   = arc(j, 2);

        arc_nwl_std = std(N_wl(arc_start:arc_end));

        for k = (arc_start + 1) : arc_end

            dmwc = 0;
            dgfc = 0; 

            dgf = L_gf(k) - L_gf(k-1);
            elv = EA.GLOel(k, i); 

            me = 1 + (10 * exp(-elv / 10));
            smg = sig0 * me;
            thresh_gf = (4 * smg) + dl;

            if abs(dgf) > thresh_gf
                dgfc = 1; 
            end

            window_start = max(arc_start, k - 30);

            if (k - window_start) >= 3
                past_nwl_data = N_wl(window_start : k-1);
                mmw = mean(past_nwl_data);
                smw = std(past_nwl_data);
                smw = max(smw, 0.01);

                dmw = N_wl(k) - mmw;

                if abs(dmw) > (5 * smw)
                    dmwc = 1; 
                end
            end

            if dgfc || (dmwc && arc_nwl_std > 0.6)

                delta_Nwl = N_wl(k) - N_wl(k-1);
                delta_Lgf = L_gf(k) - L_gf(k-1);

                A = [1, -1; lamda1, -lamda2];
                L = [delta_Nwl; delta_Lgf];
                dN_float = A \ L;

                dN1_fix = round(dN_float(1));
                dN2_fix = round(dN_float(2));

                if dN1_fix ~= 0 || dN2_fix ~= 0
                    obs.GLOL1P(k:arc_end, i) = obs.GLOL1P(k:arc_end, i) - dN1_fix;
                    obs.GLOL2P(k:arc_end, i) = obs.GLOL2P(k:arc_end, i) - dN2_fix;

                    L_gf(k:arc_end) = L_gf(k:arc_end) - (dN1_fix * lamda1 - dN2_fix * lamda2);
                    N_wl(k:arc_end) = N_wl(k:arc_end) - (dN1_fix - dN2_fix);
                end
            end
        end
    end
end
end


function [obs] = BDS_slip(obs, EA)
% INPUT:
%      obs: struct of rinex observation files, containing fields like
%           BDSC2I, BDSC6I, BDSL2I, BDSL6I.
%      EA:  struct containing elevation angles. Expected field: EA.BDSel.
%
% OUTPUT:
%      obs: Modified observation data with cycle slips repaired.

[epoch_num, sat_num] = size(obs.BDSC2I);
c = 299792458;

dt = 30; 
sig0 = sqrt(2 * (0.0027^2 + 0.0017^2)); 

iono_rate_per_hour = 4; 
dl = (iono_rate_per_hour) * (dt / 3600);

f2 = 1561.098e6; 
f6 = 1268.520e6; 

lamda2 = c / f2;
lamda6 = c / f6;
lamda_w = c / (f2 - f6); 

valid_mask = obs.BDSL2I ~= 0 & obs.BDSL6I ~= 0 & obs.BDSC2I ~= 0 & obs.BDSC6I ~= 0;

L_gf = zeros(epoch_num, sat_num);
L2_m = obs.BDSL2I .* lamda2;
L6_m = obs.BDSL6I .* lamda6;
L_gf(valid_mask) = L2_m(valid_mask) - L6_m(valid_mask);

N_wl = zeros(epoch_num, sat_num);
L_wl = (f2 * L2_m - f6 * L6_m) / (f2 - f6);
P_nl = (f2 * obs.BDSC2I + f6 * obs.BDSC6I) / (f2 + f6);

P_nl(P_nl == 0) = NaN; 
N_wl(valid_mask) = (L_wl(valid_mask) - P_nl(valid_mask)) / lamda_w;

for i = 1:sat_num 

    arc = Get_arc(valid_mask(:, i));
    if isempty(arc)
        continue;
    end

    arc( (arc(:, 2) - arc(:, 1) + 1) < 30, :) = [];
    if isempty(arc)
        continue;
    end

    for j = 1:size(arc, 1) 
        arc_start = arc(j, 1);
        arc_end   = arc(j, 2);

        arc_nwl_std = std(N_wl(arc_start:arc_end, i));

        for k = (arc_start + 1) : arc_end

            dmwc = 0; 
            dgfc = 0; 

            dgf = L_gf(k, i) - L_gf(k-1, i);
            elv = EA.BDSel(k, i); 

            me = 1 + (10 * exp(-elv / 10));
            smg = sig0 * me;
            thresh_gf = (4 * smg) + dl;

            if abs(dgf) > thresh_gf
                dgfc = 1; 
            end

            window_start = max(arc_start, k - 30);

            if (k - window_start) >= 3
                past_nwl_data = N_wl(window_start : k-1, i);
                mmw = mean(past_nwl_data);
                smw = std(past_nwl_data);
                smw = max(smw, 0.01);

                dmw = N_wl(k, i) - mmw;

                if abs(dmw) > (5 * smw)
                    dmwc = 1; 
                end
            end

            if dgfc || (dmwc && arc_nwl_std > 0.6)

                delta_Nwl = N_wl(k, i) - N_wl(k-1, i);
                delta_Lgf = L_gf(k, i) - L_gf(k-1, i);

                A = [1, -1; lamda2, -lamda6];
                L = [delta_Nwl; delta_Lgf];
                dN_float = A \ L;

                dN2_fix = round(dN_float(1)); 
                dN6_fix = round(dN_float(2)); 

                if dN2_fix ~= 0 || dN6_fix ~= 0
               
                    obs.BDSL2I(k:arc_end, i) = obs.BDSL2I(k:arc_end, i) - dN2_fix;
                    obs.BDSL6I(k:arc_end, i) = obs.BDSL6I(k:arc_end, i) - dN6_fix;

                    L_gf(k:arc_end, i) = L_gf(k:arc_end, i) - (dN2_fix * lamda2 - dN6_fix * lamda6);
                    N_wl(k:arc_end, i) = N_wl(k:arc_end, i) - (dN2_fix - dN6_fix);
                end
            end
        end
    end
end
end



function [obs] = BDS1_slip(obs, EA)
% INPUT:
%      obs: struct of rinex observation files, containing fields like
%           BDSC1I, BDSC7I, BDSL1I, BDSL7I.
%      EA:  struct containing elevation angles. Expected field: EA.BDSel.
%
% OUTPUT:
%      obs: Modified observation data with cycle slips repaired.

[epoch_num, sat_num] = size(obs.BDSC1I);
c = 299792458; 

dt = 30; 
sig0 = sqrt(2 * (0.0027^2 + 0.0017^2));

iono_rate_per_hour = 4; 
dl = (iono_rate_per_hour) * (dt / 3600);

f_B1I = 1561.098e6;
f_B3I = 1207.140e6;

lamda_B1I = c / f_B1I;
lamda_B3I = c / f_B3I;
lamda_w = c / (f_B1I - f_B3I); 

valid_mask = obs.BDSL1I ~= 0 & obs.BDSL7I ~= 0 & obs.BDSC1I ~= 0 & obs.BDSC7I ~= 0;

L_gf = zeros(epoch_num, sat_num);
L1I_m = obs.BDSL1I .* lamda_B1I;
L7I_m = obs.BDSL7I .* lamda_B3I;
L_gf(valid_mask) = L1I_m(valid_mask) - L7I_m(valid_mask);

N_wl = zeros(epoch_num, sat_num);
L_wl = (f_B1I * L1I_m - f_B3I * L7I_m) / (f_B1I - f_B3I);
P_nl = (f_B1I * obs.BDSC1I + f_B3I * obs.BDSC7I) / (f_B1I + f_B3I);

P_nl(P_nl == 0) = NaN; 
N_wl(valid_mask) = (L_wl(valid_mask) - P_nl(valid_mask)) / lamda_w;

for i = 1:sat_num 
   
    arc = Get_arc(valid_mask(:, i));
    if isempty(arc)
        continue;
    end

    arc( (arc(:, 2) - arc(:, 1) + 1) < 30, :) = [];
    if isempty(arc)
        continue;
    end

    for j = 1:size(arc, 1) 
        arc_start = arc(j, 1);
        arc_end   = arc(j, 2);

        arc_nwl_std = std(N_wl(arc_start:arc_end, i));
        
        for k = (arc_start + 1) : arc_end
            
            dmwc = 0; 
            dgfc = 0; 

            dgf = L_gf(k, i) - L_gf(k-1, i);
            elv = EA.BDSel(k, i); 
            
            me = 1 + (10 * exp(-elv / 10));
            smg = sig0 * me;
            thresh_gf = (4 * smg) + dl;
            
            if abs(dgf) > thresh_gf
                dgfc = 1;
            end

            window_start = max(arc_start, k - 30);
            
            if (k - window_start) >= 3
                past_nwl_data = N_wl(window_start : k-1, i);
                mmw = mean(past_nwl_data);
                smw = std(past_nwl_data);
                smw = max(smw, 0.01);
                
                dmw = N_wl(k, i) - mmw;
                
                if abs(dmw) > (5 * smw)
                    dmwc = 1;
                end
            end

            if dgfc || (dmwc && arc_nwl_std > 0.6)
                
                delta_Nwl = N_wl(k, i) - N_wl(k-1, i);
                delta_Lgf = L_gf(k, i) - L_gf(k-1, i);

                A = [1, -1; lamda_B1I, -lamda_B3I];
                L = [delta_Nwl; delta_Lgf];
                dN_float = A \ L;
                
                dN_B1I_fix = round(dN_float(1)); 
                dN_B3I_fix = round(dN_float(2)); 
                
                if dN_B1I_fix ~= 0 || dN_B3I_fix ~= 0
                    obs.BDSL1I(k:arc_end, i) = obs.BDSL1I(k:arc_end, i) - dN_B1I_fix;
                    obs.BDSL7I(k:arc_end, i) = obs.BDSL7I(k:arc_end, i) - dN_B3I_fix;

                    L_gf(k:arc_end, i) = L_gf(k:arc_end, i) - (dN_B1I_fix * lamda_B1I - dN_B3I_fix * lamda_B3I);
                    N_wl(k:arc_end, i) = N_wl(k:arc_end, i) - (dN_B1I_fix - dN_B3I_fix);
                end
            end
        end
    end
end
end



function [obs] = GAL_slip(obs, EA)
% INPUT:
%      obs: struct of rinex observation files, containing fields like
%           GALC1C, GALL1C, GALC5Q, GALL5Q.
%      EA:  struct containing elevation angles. Expected field: EA.GALel.
%
% OUTPUT:
%      obs: Modified observation data with cycle slips repaired.

[epoch_num, sat_num] = size(obs.GALC1C);
c = 299792458; 

dt = 30; 
sig0 = sqrt(2 * (0.0027^2 + 0.0017^2));

iono_rate_per_hour = 4; 
dl = (iono_rate_per_hour) * (dt / 3600);

f1 = 1575.42e6;  
f5 = 1176.45e6; 

lamda1 = c / f1;
lamda5 = c / f5;
lamda_w = c / (f1 - f5); 

valid_mask = obs.GALL1C ~= 0 & obs.GALL5Q ~= 0 & obs.GALC1C ~= 0 & obs.GALC5Q ~= 0;

L_gf = zeros(epoch_num, sat_num);
L1_m = obs.GALL1C .* lamda1;
L5_m = obs.GALL5Q .* lamda5;
L_gf(valid_mask) = L1_m(valid_mask) - L5_m(valid_mask);

N_wl = zeros(epoch_num, sat_num);
L_wl = (f1 * L1_m - f5 * L5_m) / (f1 - f5);
P_nl = (f1 * obs.GALC1C + f5 * obs.GALC5Q) / (f1 + f5);

P_nl(P_nl == 0) = NaN; 
N_wl(valid_mask) = (L_wl(valid_mask) - P_nl(valid_mask)) / lamda_w;

for i = 1:sat_num 

    arc = Get_arc(valid_mask(:, i));
    if isempty(arc)
        continue;
    end

    arc( (arc(:, 2) - arc(:, 1) + 1) < 30, :) = [];
    if isempty(arc)
        continue;
    end

    for j = 1:size(arc, 1) 
        arc_start = arc(j, 1);
        arc_end   = arc(j, 2);

        arc_nwl_std = std(N_wl(arc_start:arc_end, i));

        for k = (arc_start + 1) : arc_end

            dmwc = 0; 
            dgfc = 0; 

            dgf = L_gf(k, i) - L_gf(k-1, i);
            elv = EA.GALel(k, i); 

            me = 1 + (10 * exp(-elv / 10));
            smg = sig0 * me;
            thresh_gf = (4 * smg) + dl;

            if abs(dgf) > thresh_gf
                dgfc = 1;
            end

            window_start = max(arc_start, k - 30);

            if (k - window_start) >= 3
                past_nwl_data = N_wl(window_start : k-1, i);
                mmw = mean(past_nwl_data);
                smw = std(past_nwl_data);
                smw = max(smw, 0.01);

                dmw = N_wl(k, i) - mmw;

                if abs(dmw) > (5 * smw)
                    dmwc = 1;
                end
            end

            if dgfc || (dmwc && arc_nwl_std > 0.6)

                delta_Nwl = N_wl(k, i) - N_wl(k-1, i);
                delta_Lgf = L_gf(k, i) - L_gf(k-1, i);

                A = [1, -1; lamda1, -lamda5];
                L = [delta_Nwl; delta_Lgf];
                dN_float = A \ L;

                dN1_fix = round(dN_float(1)); 
                dN5_fix = round(dN_float(2)); 

                if dN1_fix ~= 0 || dN5_fix ~= 0
                    obs.GALL1C(k:arc_end, i) = obs.GALL1C(k:arc_end, i) - dN1_fix;
                    obs.GALL5Q(k:arc_end, i) = obs.GALL5Q(k:arc_end, i) - dN5_fix;

                    L_gf(k:arc_end, i) = L_gf(k:arc_end, i) - (dN1_fix * lamda1 - dN5_fix * lamda5);
                    N_wl(k:arc_end, i) = N_wl(k:arc_end, i) - (dN1_fix - dN5_fix);
                end
            end
        end
    end
end
end



function [obs] = GALX_slip(obs, EA)
% INPUT:
%      obs: struct of rinex observation files, containing fields like
%           GALC1X, GALL1X, GALC5X, GALL5X.
%      EA:  struct containing elevation angles. Expected field: EA.GALel.
%
% OUTPUT:
%      obs: Modified observation data with cycle slips repaired.

[epoch_num, sat_num] = size(obs.GALC1X);
c = 299792458;

dt = 30; 

sig0 = sqrt(2 * (0.0027^2 + 0.0017^2));

iono_rate_per_hour = 4;
dl = (iono_rate_per_hour) * (dt / 3600);

f1 = 1575.42e6;  
f5 = 1176.45e6;  

lamda1 = c / f1;
lamda5 = c / f5;
lamda_w = c / (f1 - f5); 

valid_mask = obs.GALL1X ~= 0 & obs.GALL5X ~= 0 & obs.GALC1X ~= 0 & obs.GALC5X ~= 0;

L_gf = zeros(epoch_num, sat_num);
L1_m = obs.GALL1X .* lamda1;
L5_m = obs.GALL5X .* lamda5;
L_gf(valid_mask) = L1_m(valid_mask) - L5_m(valid_mask);

N_wl = zeros(epoch_num, sat_num);
L_wl = (f1 * L1_m - f5 * L5_m) / (f1 - f5);
P_nl = (f1 * obs.GALC1X + f5 * obs.GALC5X) / (f1 + f5);
P_nl(P_nl == 0) = NaN; 
N_wl(valid_mask) = (L_wl(valid_mask) - P_nl(valid_mask)) / lamda_w;

for i = 1:sat_num 

    arc = Get_arc(valid_mask(:, i));
    if isempty(arc)
        continue;
    end

    arc( (arc(:, 2) - arc(:, 1) + 1) < 30, :) = [];
    if isempty(arc)
        continue;
    end

    for j = 1:size(arc, 1) 
        arc_start = arc(j, 1);
        arc_end   = arc(j, 2);

        arc_nwl_std = std(N_wl(arc_start:arc_end, i));

        for k = (arc_start + 1) : arc_end

            dmwc = 0; 
            dgfc = 0; 

            dgf = L_gf(k, i) - L_gf(k-1, i);
            elv = EA.GALel(k, i);

            me = 1 + (10 * exp(-elv / 10));
            smg = sig0 * me;
            thresh_gf = (4 * smg) + dl;

            if abs(dgf) > thresh_gf
                dgfc = 1;
            end

            window_start = max(arc_start, k - 30);

            if (k - window_start) >= 3
                past_nwl_data = N_wl(window_start : k-1, i);
                mmw = mean(past_nwl_data);
                smw = std(past_nwl_data);
                smw = max(smw, 0.01);

                dmw = N_wl(k, i) - mmw;

                if abs(dmw) > (5 * smw)
                    dmwc = 1;
                end
            end

            if dgfc || (dmwc && arc_nwl_std > 0.6)

                delta_Nwl = N_wl(k, i) - N_wl(k-1, i);
                delta_Lgf = L_gf(k, i) - L_gf(k-1, i);

                A = [1, -1; lamda1, -lamda5];
                L = [delta_Nwl; delta_Lgf];
                dN_float = A \ L;

                dN1_fix = round(dN_float(1));
                dN5_fix = round(dN_float(2)); 

                if dN1_fix ~= 0 || dN5_fix ~= 0

                    obs.GALL1X(k:arc_end, i) = obs.GALL1X(k:arc_end, i) - dN1_fix;
                    obs.GALL5X(k:arc_end, i) = obs.GALL5X(k:arc_end, i) - dN5_fix;

                    L_gf(k:arc_end, i) = L_gf(k:arc_end, i) - (dN1_fix * lamda1 - dN5_fix * lamda5);
                    N_wl(k:arc_end, i) = N_wl(k:arc_end, i) - (dN1_fix - dN5_fix);
                end
            end
        end
    end
end
end


%%
function arc = Get_arc(array)
len=length(array);
arc=[];
for i=1:len
    if i==len
        if array(i)~=0
            arc=[arc,i];
        end
        continue;
    end
    if i==1&&array(i)~=0
        arc=[arc,i];
    end
    if array(i)==0&&array(i+1)~=0
        arc=[arc,i+1];
        continue;
    end
    if array(i)~=0&&array(i+1)==0
        arc=[arc,i];
        continue;
    end
end
if len==0,return,end
arc=reshape(arc,2,[]);
arc=arc';
end
%%  Observation file corresponding file
function new_files=now_files(stations,files)
file_features = cell(length(files)-2, 1);  
for i = 1:length(files)
    file_name = files(i).name;
    file_features{i} = file_name(1:4);  
end
valid_files_indices = [];
for i = 1:length(file_features)
    if ismember(file_features{i}, stations)
        valid_files_indices = [valid_files_indices, i];
    end
end
new_files = files(valid_files_indices);
end
