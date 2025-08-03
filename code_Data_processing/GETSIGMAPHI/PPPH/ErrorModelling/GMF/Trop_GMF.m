function [Trop,Mwet,Mn,Me,ZHD] = Trop_GMF(rec,sat,dmjd,p)
%% Calculate Tropospheric Delay using GMF Model
% Function:
%     Calculates the tropospheric delay correction for GNSS signals using the
%     Global Mapping Function (GMF) model. This function computes the slant
%     hydrostatic delay and provides the wet and gradient mapping functions.
%
% INPUT:
%     rec:        Receiver position vector in ECEF [X; Y; Z], in meters.
%     sat:        Satellite position vector in ECEF [X; Y; Z], in meters.
%     dmjd:       Modified Julian Date.
%     p:          Atmospheric pressure at the receiver site, in hPa.
%
% OUTPUT:
%     Trop:       Slant Hydrostatic Delay (SHD), in meters.
%     Mwet:       GMF wet mapping function.
%     Mn:         North component of the tropospheric gradient mapping function.
%     Me:         East component of the tropospheric gradient mapping function.
%     ZHD:        Zenith Hydrostatic Delay, in meters.
%% ---------------------------------------------------------------------
[ellp] = xyz2plh(rec,0);
dlat = ellp(1);
dlon = ellp(2);
hell = ellp(3);

[Az,Elv] = local(rec,sat,0);

f = 0.0022768;
k = 1 - (0.00266*cos(2*dlat)) - (0.28*10^-6*hell);
ZHD = f*(p/k);

[gmfh,gmfw] = gmf_f_hu(dmjd,dlat,dlon,hell,(pi/2 - Elv));

Trop = gmfh*ZHD;
Mwet = gmfw;
Mg = 1/((tan(Elv)*sin(Elv))+0.0032);
Mn = Mg*cos(Az);
Me = Mg*sin(Az);
end

