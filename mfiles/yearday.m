function [yday]=yearday(day,month,year,hours,minutes,sec)
% Usage: [yday]=yearday(day,month,year,hours,minutes,sec), e.g. year_day(3,7,1994,7,40,3)
%  inputs:
%    day: day of month, e.g. 1 for first day, 31 for last day of August
%    month: month as an integer
%    year:
%    hours: integer hours elapsed in this day, 0 for 1230 am, 23 for 1130 pm 
%    minutes: integer minutes elapsed in this hour
%    sec: decimal seconds elapsed in this minute
%  output:
%    yday: Fraction of year elapsed to present moment.  Hence,
%     1201 am on 31 Jan is yday=1/(365*24*60) when the new year is 
%     not a leap year.
%  M.Gregg, revised 16sep96 to make  1 Jan day 0.  Previously,
%    1 Jan was day 1.

% Set up vector with the number of days in each month
month_days=[31 28 31 30 31 30 31 31 30 31 30 31];
Leap = 0;  % Use first valid year to test whether leap year
ix = find(year>1800);
if ~isempty(ix)
    ix=ix(1);
else
    ix=1;
end
if (rem(year(ix),4)==0 & year(ix) ~= 2000) | rem(year(ix),400)==0
	month_days(2)=29;
    Leap = 1;
end

% Calc days elapsed at start of each month
cumdays=cumsum(month_days); 
cumdays=[0 cumdays(1:11)]; % 

% Calc # of complete days elapsed
imok = find(month>=1 & month<=12);
yday = NaN*month;
yday(imok) = cumdays(month(imok)) + day(imok)-1;

% Add fraction of day elapsed	
yday = yday + hours/24 + minutes/(24*60) + sec/(24*3600);

