function [wup] = wind_up(rec,sat,sun,prev)
%% Calculate Phase Wind-Up Correction
% Function:
%     Calculates the carrier phase wind-up effect, which arises from the
%     relative rotation of the satellite and receiver antennas.
%
% INPUT:
%     rec:        Receiver position vector in ECEF [X; Y; Z], in meters.
%     sat:        Satellite position vector in ECEF [X; Y; Z], in meters.
%     sun:        Sun position vector in ECEF [X; Y; Z], in meters.
%     prev:       The wind-up correction from the previous epoch, in radians,
%                 used to maintain phase continuity.
%
% OUTPUT:
%     wup:        Calculated phase wind-up correction for the current epoch,
%                 in radians.
%% ---------------------------------------------------------------------
esun = sun - sat;
esun = esun./norm(esun);
ez   = -1.*sat; ez = ez./norm(ez);
ey   = cross(ez,esun); ey = ey./norm(ey);
ex   = cross(ey,ez);   ex = ex./norm(ex);
xs   = ex; ys = ey;

[elip] = xyz2plh(rec,0);
phi = elip(1);
lam = elip(2);

xr   = [(-sin(phi)*cos(lam)) (-sin(phi)*sin(lam)) cos(phi)];
yr   = [  sin(lam) -cos(lam) 0];

k = rec - sat; k = k./norm(k);

Ds = xs - k.*(dot(k,xs)) - cross(k,ys);
Dr = xr - k.*(dot(k,xr)) + cross(k,yr);
wup= acos(dot(Ds,Dr)/norm(Ds)/norm(Dr));
if dot(k,(cross(Ds,Dr)))<0, wup = -wup; end
wup = (2*pi*floor(((prev - wup)/2/pi)+0.5)) + wup;
end

