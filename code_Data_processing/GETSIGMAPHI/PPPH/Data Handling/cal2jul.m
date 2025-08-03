function [jd,mjd] = cal2jul(year,mon,day,sec)
%% Convert Calendar Date to Julian Date
% Function:
%     Converts a Gregorian calendar date (year, month, day) and time (as
%     seconds of the day) to Julian Date (JD) and Modified Julian Date (MJD).
%     The Julian Date is a continuous count of days and is widely used in
%     astronomy, geodesy, and other scientific applications for timekeeping.
%
% INPUT:
%     year:       Year (scalar).
%     mon:        Month of the year (1-12, scalar).
%     day:        Day of the month (1-31, scalar).
%     sec:        Seconds of the day (elapsed seconds since 00:00 of the given date).
%
% OUTPUT:
%     jd:         Julian Date.
%     mjd:        (Optional) Modified Julian Date (MJD = JD - 2400000.5).
%
% Note:
%     The function includes basic checks for the number and range of input arguments.
%% ---------------------------------------------------------------------
narginchk(4,4)

sec = sec/3600;

if ~isscalar(year) || ~isscalar(mon) || ~isscalar(day)
    error('Year, Month and Day should be scalar.')
end

if mon<1 || mon>12
    error('Month should be between 1 and 12.')
end


if day<1 || day>31
    error('Day should be between 1 and 31.')
end

if mon<=2
    m = mon + 12;
    y = year - 1;
else
    m = mon;
    y = year;
end

jd = floor(365.25*y) + floor(30.6001*(m+1)) + day + (sec/24) + 1720981.5;
mjd= jd - 2400000.5;

nargoutchk(1,2)
end

