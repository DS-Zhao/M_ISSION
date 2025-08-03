function [rapc] = rec_apc(s_xyz,r_xyz,r_apc,k)
%% Calculate Receiver Antenna Phase Center (APC) Correction
% Function:
%     Calculates the receiver Antenna Phase Center (APC) correction by
%     projecting the APC offset onto the line-of-sight (LOS) vector.
%
% INPUT:
%     s_xyz:      Position vector of the source (e.g., satellite) in ECEF
%                 [X; Y; Z], in meters.
%     r_xyz:      Position vector of the receiver in ECEF [X; Y; Z], in meters.
%     r_apc:      A 3D array containing receiver APC offsets in the local
%                 ENU (East, North, Up) frame, in meters.
%     k:          Index to select the appropriate APC offset from r_apc.
%
% OUTPUT:
%     rapc:       Scalar receiver APC correction value, in meters.
%% ---------------------------------------------------------------------
l   = r_xyz - s_xyz;
los = l./norm(l);

[elip] = xyz2plh(r_xyz,0);
lat = elip(1);
lon = elip(2);

ori = [(-sin(lon)) (-cos(lon)*sin(lat)) (cos(lon)*cos(lat));...
        (cos(lon)) (-sin(lon)*sin(lat)) (sin(lon)*cos(lat));...
               (0) (cos(lat)) (sin(lat))];
           
f = r_apc(:,:,k);

c = [f(2);
     f(1);
     f(3)];

p = ori*c;
if size(los,1)~=size(p,1)
    los = los';
end

rapc = dot(p,los);

end

