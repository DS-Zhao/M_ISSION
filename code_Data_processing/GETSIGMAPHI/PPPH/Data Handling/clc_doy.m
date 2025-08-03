function [doy] = clc_doy(year,mon,day)
%% Calculate Day of Year (DOY)
% Function:
%     Calculates the sequential day number within a given year (Day of Year, DOY)
%     from a standard calendar date (year, month, day).
%     The function first determines whether the input year is a leap year to
%     correctly establish the number of days in February. It then computes the
%     DOY by summing the days of all preceding months and adding the day of the
%     current month.
%
% INPUT:
%     year:       The year.
%     mon:        The month of the year (1-12).
%     day:        The day of the month (1-31).
%
% OUTPUT:
%     doy:        The Day of Year (ranging from 1 to 365 for a common year,
%                 or 1 to 366 for a leap year).
%
% Example:
%     clc_doy(2024, 2, 10) will return 41 (since 2024 is a leap year, 31+10=41)
%% ---------------------------------------------------------------------


if rem(year,400)==0
    sit = 1;
elseif rem(year,300)==0
    sit = 0;
elseif rem(year,200)==0
    sit = 0;
elseif rem(year,100)==0
    sit = 0;
elseif rem(year,4)==0
    sit = 1;
else
    sit = 0;
end

switch sit
    case 0
        dom = [31 28 31 30 31 30 31 31 30 31 30 31];
    case 1
        dom = [31 29 31 30 31 30 31 31 30 31 30 31];
end

if mon-1<1
    doy = day;
else
    doy = sum(dom(1:(mon-1))) + day;
end

end

