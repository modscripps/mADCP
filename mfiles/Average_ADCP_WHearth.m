function Vel = Average_ADCP_WHearth(AllVel, AvgInt, TimeVecInt)

% Average_ADCP_WHearth Average moored ADCP data in time
%
%   Vel = Average_ADCP_WHearth(AllVel, AvgInt, TimeVecInt)
%
%   INPUT   AllVel -     Single ping ADCP data in earth coordinates (output
%                        from Run_WHbeams_ALL)
%           AvgInt -     Averaging interval (in hours)
%           TimeVecInt - Interval of new time grid (in hours)
%       
%   OUTPUT  Vel - Structure with time-averaged velocities
%
%   Gunnar Voet
%   gvoet@ucsd.edu
%
%   Changed: 2017-05-04

Vel = AllVel;

% Create new time vector
TMEgrid  = floor(AllVel.dtnum(1)*24)/24:TimeVecInt/1440:ceil(AllVel.dtnum(end)*24)/24;

SampInt = nanmean(diff(AllVel.dtnum)).*86400; % sampling interval in seconds

NNB = ceil((AvgInt.*60)./SampInt); %rounding up...this is the number of samples in the average.

fprintf('\n')
fprintf('Applying boxcar running mean\n')
uAve   = g_boxcar_smoothNONAN(AllVel.u,NNB);
vAve   = g_boxcar_smoothNONAN(AllVel.v,NNB);
wAve   = g_boxcar_smoothNONAN(AllVel.w,NNB);
ecMAve = g_boxcar_smoothNONAN(AllVel.ecAve,NNB);
tempAve = g_boxcar_smoothNONAN(AllVel.temp,NNB);
xducerdepthAve = g_boxcar_smoothNONAN(AllVel.xducer_depth,NNB);

% smoothed this earlier (but at different time scale!)
% xducer_depth = AllVel.depth_xducer;


%NaNing beginning and ends where averaging causes problems. Being overly
%conservative here.
uAve(:,1:NNB)       = NaN;
uAve(:,end-NNB:end) = NaN;
vAve(:,1:NNB)       = NaN;
vAve(:,end-NNB:end) = NaN;
wAve(:,1:NNB)       = NaN;
wAve(:,end-NNB:end) = NaN;
ecMAve(:,1:NNB)     = NaN;
ecMAve(:,end-NNB:end) = NaN;

xducerdepthAve(end-NNB:end) = NaN;
xducerdepthAve(1:NNB)       = NaN;
tempAve(end-NNB:end)      = NaN;
tempAve(1:NNB)            = NaN;
DnumAve                   = AllVel.dtnum;



%% Interpolate to new time grid

Vel.u = nan(length(Vel.z),length(TMEgrid));
Vel.v = Vel.u;
Vel.w = Vel.u;
Vel.ecAve = Vel.u;

fprintf('Interpolating to new time grid\n')
for ii=1:length(Vel.z);
Vel.u(ii,:) = interp1(DnumAve,uAve(ii,:),TMEgrid,'linear');
Vel.v(ii,:) = interp1(DnumAve,vAve(ii,:),TMEgrid,'linear');
Vel.w(ii,:) = interp1(DnumAve,wAve(ii,:),TMEgrid,'linear');
Vel.ecAve(ii,:) = interp1(DnumAve,ecMAve(ii,:),TMEgrid,'linear');
end
Vel.temp = interp1(DnumAve,tempAve,TMEgrid,'linear');
Vel.xducer_depth = interp1(DnumAve,xducerdepthAve,TMEgrid,'linear');
Vel.yday = interp1(DnumAve,Vel.yday,TMEgrid,'linear');
Vel.dtnum = TMEgrid;
fprintf('\n')
























% 
% 
% 
% 
%     
%     
%     % Creating matrices
%     if ~exist('uMat','var')==1 && ~isempty(ADP);
%         uMat  = NaN.*ones(length(ADP.depth),length(TMEgrid));
%         vMat  = uMat;
%         wMat  = uMat;
%         ecMat = uMat;
%         xDucerDep = NaN.*ones(1,length(TMEgrid));
%         tempM = xDucerDep;
%     end
% 
%     %%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%
%     %%% must have diff't approaches for ours and Jody's ADCPs..
% 
%     if InP.AVG.flag=='Y';
% 
%     SampInt = nanmean(diff(ADP.dtnum)).*86400; % sampling interval in seconds
%     AVGINT = InP.AVG.int; % averaging time in minutes.
% 
%     NNB = ceil((AVGINT.*60)./SampInt); %rounding up...this is the number of samples in the average.
% 
%     %     BB=boxcar(NNB)./nansum(boxcar(NNB));
% 
%     uAve   = g_boxcar_smoothNONAN(ADP.u,NNB);
%     vAve   = g_boxcar_smoothNONAN(ADP.v,NNB);
%     wAve   = g_boxcar_smoothNONAN(ADP.w,NNB);
%     ecMAve = g_boxcar_smoothNONAN(ecM,NNB);
% 
%     %smoothed this earlier
%     xducer_depth = ADP.depth_xducer;
%     tempAve = boxcar_smoothNONAN(ADP.temp,NNB);
% 
%     %NaNing beginning and ends where averaging causes problems. Being overly
%     %conservative here.
%     uAve(:,1:NNB)       = NaN;
%     uAve(:,end-NNB:end) = NaN;
%     vAve(:,1:NNB)       = NaN;
%     vAve(:,end-NNB:end) = NaN;
%     wAve(:,1:NNB)       = NaN;
%     wAve(:,end-NNB:end) = NaN;
%     ecMAve(:,1:NNB)     = NaN;
%     ecMAve(:,end-NNB:end) = NaN;
% 
%     xducer_depth(end-NNB:end) = NaN;
%     xducer_depth(1:NNB)       = NaN;
%     tempAve(end-NNB:end)      = NaN;
%     tempAve(1:NNB)            = NaN;
%     DnumAve                   = ADP.dtnum;
% 
%     else % no averaging
% 
%     uAve=ADP.u;
%     vAve=ADP.v;
%     wAve=ADP.w;
%     ecMAve=ecM;
%     xducer_depth=ADP.depth_xducer;
%     tempAve=ADP.temp;
%     DnumAve=ADP.dtnum;
% 
%     end
% 
%     % To save time only interpolating onto a portion of the total time.
%     Isff = find(TMEgrid>=DnumAve(1) & TMEgrid<=DnumAve(end)); 
% 
%     % Interpolating onto standard grid.
%     [~,Iuu] = unique(DnumAve);
% 
%     for ii=1:length(ADP.depth);
% 
%       %u
%       uAveI = interp1(DnumAve(Iuu),uAve(ii,Iuu),TMEgrid(Isff),'linear'); 
%       iso   = find(~isnan(uAveI));
%       uMat(ii,Isff(iso)) = uAveI(iso);
% 
%       %v
%       vAveI = interp1(DnumAve(Iuu),vAve(ii,Iuu),TMEgrid(Isff),'linear');
%       iso   = find(~isnan(vAveI));
%       vMat(ii,Isff(iso)) = vAveI(iso);
% 
%       %w
%       wAveI = interp1(DnumAve(Iuu),wAve(ii,Iuu),TMEgrid(Isff),'linear'); 
%       iso   = find(~isnan(wAveI));
%       wMat(ii,Isff(iso)) = wAveI(iso);
% 
%       %ec
%       ecAveI = interp1(DnumAve(Iuu),ecMAve(ii,Iuu),TMEgrid(Isff),'linear'); 
%       iso    = find(~isnan(ecAveI));
%       ecMat(ii,Isff(iso)) = ecAveI(iso);
% 
%     end
% 
%     xducerDepthI = interp1(DnumAve(Iuu),xducer_depth(Iuu),TMEgrid(Isff),'linear'); 
%     iso          = find(~isnan(xducerDepthI));
%     xDucerDep(Isff(iso)) = xducerDepthI(iso);
% 
% 
%     tempAveI = interp1(DnumAve(Iuu),tempAve(Iuu),TMEgrid(Isff),'linear'); 
%     iso = find(~isnan(tempAveI));
%     tempM(Isff(iso)) = tempAveI(iso);
%   end
% % end
% 
% % finding all nans at beginning or end of timeseries, removing from
% % timeseries
% 
% %this gets all nans, just removing beginning and end sections
% Igood=find(~isnan(nanmean(uMat)));
% 
% if ~isempty(Igood);
%     IgdS=Igood(1);
%     IgdE=Igood(end);
% else
%     disp('No good data within range specified');
%     return
% end
% 
% % now applying time offset over entire timeseries..
% % TmeOff
% 
% TmeOffL = linspace(0,TmeOff./86400,length(TMEgrid)); %linear offset fix.
% TMEgridn = TMEgrid+TmeOffL;
% 
% % Assign data to output structure
% Vel.u = uMat(:,IgdS:IgdE);
% Vel.v = vMat(:,IgdS:IgdE);
% Vel.w = wMat(:,IgdS:IgdE);
% Vel.ecAve = ecMat(:,IgdS:IgdE);
% Vel.dtnum = TMEgridn(IgdS:IgdE);
% 
% ydaya = datenum2yday(TMEgridn(IgdS));
% Vel.yday = ydaya+Vel.dtnum-Vel.dtnum(1);
% 
% Vel.z = depthvec;
% Vel.temp = tempM(IgdS:IgdE);
% Vel.xducer_depth = xDucerDep(IgdS:IgdE);
% Vel.botdepth = InP.NomBotDepth;
% Vel.lat = Lat;
% Vel.lon = Lon;
% 
% Vel.MagDecUsed = MagDec;
% 
% Vel.info.snADCP = WHtyp(1:end);
% Vel.info.MooringID = MooringID;
% Vel.info.Cruise = CrName;
% Vel.info.startyear = StartYear;
% Vel.info.timeOffset = TmeOff;
% Vel.info.timeOffset_notes = 'applied. Seconds that ADCP was ahead of GMT--linearly applied over timeseries';
% 
% Vel.info.sysconfig = sysconfig;
% Vel.info.NomInstDepth = InP.NomInstDepth;
% Vel.info.mfiles=['in order: Conv_ADCP_mooring_[experiment] (coverts to mat files)'...
%     'Run_WHbeams_[experiment]_AVG.m (conv. to earth coord +avgs. & calls get_ADCP_Any_dnum.m, WHbeam_Process.m, etc.)'];
% Vel.info.processing=['10 1-s water pings & 10 BT pings avg. into ' num2str(TgridInt) ' minute ensembles.'];
% if exist('Notes','var')==1;
%     Vel.info.notes=Notes;
% end
% 
% 
% 
% 
