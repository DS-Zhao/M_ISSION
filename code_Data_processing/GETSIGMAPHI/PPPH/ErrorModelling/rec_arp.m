function [rarp] = rec_arp(s_xyz,r_xyz,arp)
%% Calculate Receiver Antenna Phase Center Offset (PCO) Correction
% Function:
%     Calculates the range correction due to the receiver's antenna phase
%     center offset (PCO). It transforms the PCO vector from the local
%     topocentric frame (ENU) to the geocentric frame (ECEF) and then
%     projects it onto the satellite-to-receiver line-of-sight vector.
%
% INPUT:
%     s_xyz: Satellite position vector in ECEF [X; Y; Z].
%     r_xyz: Receiver position vector in ECEF [X; Y; Z].
%     arp:   Receiver antenna PCO vector, ordered as [Up; East; North].
%
% OUTPUT:
%     rarp:  Range correction due to the receiver antenna PCO.
%% ---------------------------------------------------------------------

l   = r_xyz -s_xyz;
los = l./norm(l);

[elip] = xyz2plh(r_xyz,0);
lat = elip(1);
lon = elip(2);

ori = [(-sin(lon)) (-cos(lon)*sin(lat)) (cos(lon)*cos(lat));...
        (cos(lon)) (-sin(lon)*sin(lat)) (sin(lon)*cos(lat));...
               (0)           (cos(lat))          (sin(lat))];
           
enu = [arp(2);arp(3);arp(1)];
p = ori*enu;
if size(p,1)~=size(los,1)
    los = los';
end
rarp = dot(p,los);
end

