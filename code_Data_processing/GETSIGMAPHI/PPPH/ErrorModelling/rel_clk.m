function [rclk] = rel_clk(s_xyz,v_xyz)
%% Calculate Relativistic Clock Correction
% Function:
%     Calculates the relativistic clock correction due to the satellite's
%     velocity and position (Sagnac effect).
%
% INPUT:
%     s_xyz:      Satellite position vector in ECEF [X; Y; Z], in meters.
%     v_xyz:      Satellite velocity vector in ECEF [X; Y; Z], in meters/second.
%
% OUTPUT:
%     rclk:       The relativistic clock correction converted to a range error,
%                 in meters.
%% ---------------------------------------------------------------------
c = 299792458; %m/s

rclk = ((-2)*((dot(s_xyz,v_xyz))/(c^2)))*c;

end

