function [rpath] = rpath(r_xyz,s_xyz)
%% Calculate Relativistic Path Delay (Shapiro Delay)
% Function:
%     Calculates the relativistic correction to the signal path length due to
%     the gravitational field of the Earth (Shapiro effect), based on receiver
%     and satellite positions in a geocentric frame.
%
% INPUT:
%     r_xyz:      Receiver position vector in a geocentric frame (e.g., ECEF)
%                 [X; Y; Z], in meters.
%     s_xyz:      Satellite position vector in a geocentric frame (e.g., ECEF)
%                 [X; Y; Z], in meters.
%
% OUTPUT:
%     rpath:      Relativistic path delay, in meters.
%% ---------------------------------------------------------------------
% WGS84 constants
mu=3986004.418*(10^8);  %m3/s2
c =2.99792458*(10^8);   %m/s

rsat = norm(s_xyz);
rrec = norm(r_xyz);
rs   = s_xyz - r_xyz;
rrs  = norm(rs);
rpath = ((2*mu)/(c^2))*log((rsat+rrec+rrs)/(rsat+rrec-rrs));

end

