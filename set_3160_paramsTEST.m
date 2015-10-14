% Processing parameters for moored ADCP
% Parameters shown here are slightly modified from
% SN3160 on ArcticMix15 mooring AM1

%% CRUISE SPECIFIC

% Cruise short name, e.g. MC09, PAPA 08, this should be the same as the directory name
InP.Cruise.Name = 'TEST';
% Directory to cruise structure on computer
InP.Cruise.Dir = 'sample_data/';
% Within this directory the raw ADCP files should be in 'ADCP/SNxxx/data_raw

% Start before ADCP logging, end after mooring out of water.  %time in UTC
InP.Cruise.StartTime = '05-Sep-2015 01:00:00';
InP.Cruise.EndTime   = '05-Sep-2015 04:00:00';
InP.Cruise.StartYear = 2015;


%% MOORING SPECIFIC

% Enter mooring name here. This is used for info only
InP.MooringID = 'AM1';
InP.Lat =     72+35.646./60;
InP.Lon =  -(145+01.002./60);  % Negative for west longitude

% Get this at http://www.ngdc.noaa.gov/geomag-web/#declination
% or from a chart.
% positive is EAST declination.
InP.MagDec = 20;  % East decliniation is positive!
                  % 20degE is a good value for the arctic mooring

% Time offset in seconds at end of deployment, positive is ADCP ahead of GMT.
InP.TimeOffset = -9;
% This is a constant offset for the entire duration of the deployment
% in HOURS (positive is HOURS BEHIND GMT).
% Usually this does not need to be different from 0.
InP.TimeShift  = 0;

% Nominal bottom depth at the mooring
InP.NomBotDepth  = 3460;
% Nominal instrument depth
InP.NomInstDepth = 50;

%% INSTRUMENT SPECIFIC

% Instrument serial number.
InP.snADCP = '3160';
InP.deployno = NaN; % NaN for single deployments, 1 or more for multiple
                    % deployments of the same instrument (even if on
                    % different moorings).

% Depth grid for final product
InP.ZGrid = [0:1:50]';  %in true depth coordinates
                    
%InP.set.theta_o=20; %Angle of the ADCP beams from vertical
%InP.set.Cnvx = 1; % convex config, 0 is concave
InP.set.UpDown = 0; % down facing = 1, else up facing  %this will be checked vs. ADCP header info.

InP.set.Psrc = 1; % pressure source: 1=ADCP internal, OR 2xN vector=[yday; dbar] for interp1
InP.set.Ztyp = 1; % Type of output depth grid: surface rel = 1; ADCP rel = 0; 
InP.set.Ptyp = 0; % beam coords and all pings = 0, earth coords (ensembles) = 1;
InP.set.Sal0 = -1; % <0=use recorded soundspeed, >0(eg 35)= use this salinity for soundspeed calc
   


%% DATA CULLING / QUALITY CONTROL SETTINGS
                     
InP.qc.Xblank = 7;  % Blanks out Xblank percent of the instrument depth (uplooking) or 
                    % height above the bottom (downlooking) near the surface and bottom respectively. 
                    % default is 10 (percent).

InP.qc.WC_val = 10;  %this is the Low Correlation threshold to be used for data in beam coords. (RDI default is 64).
InP.qc.werrMax=1; %culls data above this error threshold (default is 1);
InP.qc.stdMax=4; %culls data greater than X standard deviations (2 is good) 

InP.qc.velMax=2.00; %culls data greater than x 
InP.qc.diffMax=1;  %culls data in adjacent bins (in time) with abs differences greater than this (does forward
%and back to get spikes.

InP.qc.ExcludeBins = [];

InP.AVG.flag='Y';
InP.AVG.int=5; %averaging interval in minutes when AVG.flag='Y';
InP.AVG.TimeGridInt=5; %in minutes.

%%%BELOW FLAGS for PRESSURE-averaging.  If pressure recorded every ping, it
%%%is noisy, this needs to be knocked down to properly map bins if using
%%%ADCP pressure as reference. 

InP.PAVG.flag = 'Y';
InP.PAVG.int  = 30; %averaging interval in minutes when AVG.flag='Y';


InP.mfiles = {sprintf('set_%s_params%s.m',InP.snADCP,InP.Cruise.Name);...
              'cruiseStruct_init_ADCP';...
              'Conv_ADCP_mooring.m';...
              'Get_ADCP_fullSC.m';...
              'raw2mat_SCv3_BT.m'};  %list of mfiles used

SaveName = sprintf('SN%s_%s_params_%s.mat',...
                   InP.snADCP,InP.MooringID,InP.Cruise.Name);
save(SaveName,'InP')


% TODO: InP.ZGrid % this should be calculated rather than set...using the
% bin size, ADCP depth, etc.