function [leap_sec] = leapsec(mjd)
%% Get the Number of Leap Seconds for a Given Date
% Function:
%     This function returns the cumulative number of leap seconds for a given
%     Modified Julian Date (MJD). Leap seconds represent the offset between
%     Coordinated Universal Time (UTC) and International Atomic Time (TAI), where
%     TAI - UTC = (leap seconds) + 10s. The function uses a hard-coded table
%     of leap second insertions covering the period from January 1, 2000, onwards.
%
% INPUT:
%     mjd:        A numeric value representing the Modified Julian Date.
%
% OUTPUT:
%     leap_sec:   The cumulative number of leap seconds on the given MJD.
%% ---------------------------------------------------------------------


if mjd<51544
    error('The given date have to be greater than December 1st,2000');
end

if     mjd>=51544 && mjd<53736 %01.01.2000-01.01.2006
    leap_sec = 32;
elseif mjd>=53736 && mjd<54832 %01.01.2006-01.01.2009
    leap_sec = 33;
elseif mjd>=54832 && mjd<56109 %01.01.2009-01.07.2012
    leap_sec = 34;
elseif mjd>=56109 && mjd<57204 %01.07.2012-01.07.2015
    leap_sec = 35;
elseif mjd>=57204 && mjd<57754 %01.07.2015-01.01.2017
    leap_sec = 36;
else
    leap_sec = 37;
end
    
end

