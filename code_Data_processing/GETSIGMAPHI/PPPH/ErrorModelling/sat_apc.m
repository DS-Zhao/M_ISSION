function [sapc] = sat_apc(s_xyz,r_xyz,sun_xyz,s_apc,opt,sno)
%% Calculate Satellite Antenna Phase Center (APC) Correction
% Function:
%     Calculates the satellite Antenna Phase Center (APC) correction for a
%     given signal.
%
% INPUT:
%     s_xyz:      Satellite position vector in ECEF [X; Y; Z], in meters.
%     r_xyz:      Receiver position vector in ECEF [X; Y; Z], in meters.
%     sun_xyz:    Sun position vector in ECEF [X; Y; Z], in meters.
%     s_apc:      A 3D array containing satellite APC offsets for different
%                 frequencies/signals, in meters.
%     opt:        A selector (e.g., 1 or 2) to specify which APC offset to use.
%     sno:        Satellite number (e.g., PRN), used for handling specific cases.
%
% OUTPUT:
%     sapc:       Scalar satellite APC correction, in meters.
%% ---------------------------------------------------------------------
l   = r_xyz - s_xyz;
los = l./norm(l);
% satellite -fixed coordinate system
k = (-1).*(s_xyz./(norm(s_xyz)));
rs  = sun_xyz - s_xyz;
e   = rs./(norm(rs));
j   = cross(k,e)./norm(cross(k,e));
i   = cross(j,k);
sf  = [i; j; k];

de1 = (s_apc(:,:,1))';
de2 = (s_apc(:,:,2))';

if (any(isnan(de1)) || any(isnan(de2))) && (sno>58 && sno<89) % GALILEO CONVENTIONAL
    de1 = [0.2;0;0.6];
    de2 = [0.2;0;0.6];
elseif (any(isnan(de1)) || any(isnan(de2))) && (sno>88 && sno<93) % BEIDOU CONVENTIONAL
    de1 = [0.6;0;1.1];
    de2 = [0.6;0;1.1];
end
    
if opt == 1
    rk = sf\de1;
elseif opt == 2
    rk = sf\de2;
end

sapc = dot(rk,los);
end