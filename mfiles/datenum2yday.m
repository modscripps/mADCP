function [yday,year]=datenum2yday(dnum)
%function [yday,year]=datenum2yday(dnum)
%Simple converter from datenum to yday.  Beware crossing year boundaries.
%1/15/04 MHA
%

%First figure out what year it is
ayear=1900;
toosoon=1;
%Get the datenum of the first day of that year
while toosoon==1
    ayear=ayear+1;
    dn1=datenum(ayear,1,1);
    if dnum(1) < dn1
        toosoon=0;
    end
end
year=ayear-1;

dn1=datenum(year,1,1);

yday=dnum - dn1;