function Vel = Run_WHbeams_ALL(dtnum_beg, dtnum_end, WHtyp,MooringID,Project)

% Vel = Run_WHbeams_ALL(dtnum_beg, dtnum_end, WHtyp, MooringID, Project);
%
%   Converts ping beam data to velocities in geomagmetic frame,
%       returned in structure ADP .
%
%   Inputs:
%      dtnum_beg, dtnum_end - yearday range to process, empty vectors for
%                             the whole time series
%      WHtyp     - root name of matlab ADCP type (e.g. 'SN4021'));
%      MooringID - mooring name (e.g. 'A1')
%      Project   - project name (e.g. 'ARCTIC')

% GV: using the functional form of WHbeam_Process, partly vectorized to
% increase speeds.
% GV: Use this function to read all data, only average pressure. Average
% velocities in a second step.


%% Load mat index
load([WHtyp '_' MooringID '_params_',Project,'.mat']); % Load settings for this
                                                 % instrument (InP)
% settings have the following structure:
% InP = 
%          Cruise: [1x1 struct]
%        MooringID: 'A1'
%              Lat: -44.4902
%              Lon: 152.9543
%           MagDec: 15.2000
%       TimeOffset: 2
%        TimeShift: 0
%      NomBotDepth: 4768
%     NomInstDepth: 61
%           snADCP: '22476'
%         deployno: NaN
%            ZGrid: [51x1 double]
%              set: [1x1 struct]
%               qc: [1x1 struct]
%              AVG: [1x1 struct]
%             PAVG: [1x1 struct]
%           mfiles: {5x1 cell}

%% Assign parameters
StartYear = InP.Cruise.StartYear;
MagDec    = InP.MagDec;

TmeOffT   = InP.TimeOffset; % time offset in seconds at end of deployment
                            % positive is ADCP ahead of GMT.
                            % this is OVER the entire timeseries!!
                            % above offset will be scaled by the size of
                            % the subset.
TgridInt  = InP.AVG.TimeGridInt;
WC_val    = InP.qc.WC_val;
MooringID = InP.MooringID;
Lat       = InP.Lat;
Lon       = InP.Lon;
NomDepth  = InP.NomBotDepth;
if isfield(InP,'BadBeam')
  BadBeam   = InP.BadBeam;
else
  BadBeam = [];
end

Vel = [];
depthvec = [];

% This is the serial number of the ADCP: SNxxxx
% Mnms = WHtyp;
SN = WHtyp;

% For Alford mooring group, organized in ADCP s/n subfolders ...
CrName = InP.Cruise.Name;
ADtop = [InP.Cruise.Dir '/ADCP/']; % everything is/will be under this.

% Now, define default matlab folder,index ...
MatFld  = fullfile(ADtop, SN,'data_mat'); % matlab data files folder
IndFld  = MatFld; % matlab data index folder
% MatRoot = SN; % for root names of matlab files and index
% Load the frigging index file
MatIndx = fullfile(IndFld, [SN '_' CrName '_matfiles.mat']);
load(MatIndx);

theta_o = 20; % Angle of the ADCP beams from vertical
Cnvx    = 1; % convex config

UpDown = InP.set.UpDown; % down facing = 1, else up facing

% Pr/dbar at ADCP:
Psrc = InP.set.Psrc; % 1=ADCP internal,
                     % 2xN vector=[dtnum; dbar] for interp1
                     % or string pointing to ascii file with Nx2 vector
Ztyp = InP.set.Ztyp; % Type of output depth grid:
                     % surface rel = 1; ADCP rel = 0;
Ptyp = InP.set.Ptyp; % beam coords = 0, earth coords=1;

ZGrid = InP.ZGrid;
Sal0  = InP.set.Sal0; % <0=use recorded soundspeed,
                      % >0(eg 35)= use for soundspeed calc
if isfield(InP.qc,'ExcludeBins')
  ExcludeBins = InP.qc.ExcludeBins;
else
  ExcludeBins = [];
end

% Process the whole time series if dtnum_beg and dtnum_end are left empty:
if isempty(dtnum_beg)
  dtnum_beg = min(Index.dtnum_beg);
  dtnum_end = max(Index.dtnum_end);
end


%% Scale the time offset to the actual length of the ADCP time series.

TmeOff = TmeOffT.*(dtnum_end-dtnum_beg)./(max(Index.dtnum_end)-min(Index.dtnum_beg)); %scaling time offset by subset of time-series.


%% Read raw ADCP data
fprintf('\n')
fprintf('- Reading whole ADCP time series in beam coordinates...\n')
ADRaw = get_ADCP_Any_dnum(dtnum_beg, dtnum_end, MatIndx, MatFld, 1);
sysconfig = ADRaw.sysconfig;

% Warn if no data available
if isempty(ADRaw)
  error(['Run_WHbeams_ALL: Cannot get ' SN ' data, exit'])
end

%% Assign ADCP pressures:
if ischar(Psrc)
  load(Psrc);
  if size(a,1)>size(a,2)
    a = a';
  end
  Psrc = a;
end

if size(Psrc,1)==2 && size(Psrc,1)>1 % interpolation array
  ADRaw.depth_xducer = [];
  id = find(ADRaw.dtnum>=Psrc(1,1) & ADRaw.dtnum<=Psrc(1,end));
  if ~isempty(id)
    ADRaw.depth_xducer(id) = interp1(Psrc(1,:),Psrc(2,:),ADRaw.dtnum(id));
  end
end

%% Averaging depth_xducer over our averaging interval.
if InP.PAVG.flag == 'Y';
SampInta = nanmean(diff(ADRaw.dtnum)).*86400; % sampling interval
                                              % in seconds
AVGINT = InP.PAVG.int; % averaging time in minutes.
NNBa = ceil((AVGINT.*60)./SampInta); % rounding up...this is the number
                                     % of samples in the average.

depthS = boxcar_smoothNONAN(ADRaw.depth_xducer,NNBa);
igg = find(~isnan(depthS));
[~,iTy] = unique(ADRaw.dtnum(igg));
depthSI = interp1(ADRaw.dtnum(igg(iTy)),depthS(igg(iTy)),ADRaw.dtnum,'nearest','extrap');
ADRaw.depth_xducer = depthSI;
end

%% Sound velocity
ADRaw.svel_calc = [];
if Sal0 > 0 % Compute actual soundspeeds
    s0 = ones(size(ADRaw.temp)) * Sal0;
    % comment out improper one (sal,pres in CU,MPa for MCG's)
    %ADRaw.svel_calc = sw_svel(s0/1000, ADRaw.temp, ADRaw.depth_xducer/100);
    ADRaw.svel_calc = sw_svel(s0, ADRaw.temp, ADRaw.depth_xducer);
    clear s0
end

%% Convert to earth coordinates
fprintf('- Converting to earth coordinates...\n')
ADP = WHbeam_ProcessFunctionVectorized(ADRaw,UpDown,ZGrid,theta_o,Cnvx,MagDec,WC_val,ExcludeBins,BadBeam);

% Averaging echo intensities
ecMM = ones(length(ADP.depth),length(ADP.dtnum),4);
ecMM(:,:,1) = ADP.ec1;
ecMM(:,:,2) = ADP.ec2;
ecMM(:,:,3) = ADP.ec3;
ecMM(:,:,4) = ADP.ec4;
ecM = nanmean(ecMM,3);

%% getting rid of data w/in 10% of sampling depth or near bottom
% taking 10% of depth xducer, removing this from top
if UpDown==0; %upward-looking   OK, here it is clear that z_adcp is positive for down, negative for up-looking.
  Xblank = (InP.qc.Xblank./100).*(ADP.depth_xducer);
  Ido = find(ADP.depth*ones(1,length(Xblank))<=ones(length(ADP.depth),1)*Xblank);

  ADP.u(Ido)=NaN;
  ADP.v(Ido)=NaN;
  ADP.w(Ido)=NaN;

else %downward looking
  XPer = InP.qc.Xblank;

  % Cutting off velocities below bottom and those within X% of bottom
  DiffD = NomDepth-nanmedian(ADP.depth_xducer);

  BADD = ((XPer./100)).*DiffD;
  Ido  = find(ADP.depth>(NomDepth-BADD)); 

  ADP.u(Ido,:) = NaN;
  ADP.v(Ido,:) = NaN;
  ADP.w(Ido,:) = NaN;
end

%% Depth Vector
if ~isempty(ADRaw) && isempty(depthvec);
  depthvec = ADP.depth;
end

%% Find all nans at beginning or end of timeseries and remove them
Igood = find(~isnan(nanmean(ADP.u,1)));

if ~isempty(Igood);
    IgdS=Igood(1);
    IgdE=Igood(end);
else
    disp('No good data within range specified');
    return
end

%% Find depth levels without data and remove them as well
% Keep +/-1 level for nicer data display
Igood2 = find(~isnan(nanmean(ADP.u,2)));
if ~isempty(Igood2);
    IgdS2=Igood2(1);
    if IgdS2>1
      IgdS2 = IgdS2-1;
    end
    IgdE2=Igood2(end);
    if IgdE2<length(depthvec)
      IgdE2 = IgdE2+1;
    end
else
    disp('No good data within range specified');
    return
end


%% Apply time offset over entire timeseries
TmeOffL = linspace(0,TmeOff./86400,length(ADRaw.dtnum)); %linear offset fix.
TMEgridn = ADRaw.dtnum+TmeOffL;

%% Assign data to output structure
Vel.u = ADP.u(IgdS2:IgdE2,IgdS:IgdE);
Vel.v = ADP.v(IgdS2:IgdE2,IgdS:IgdE);
Vel.w = ADP.w(IgdS2:IgdE2,IgdS:IgdE);

Vel.ecAve = ecM(IgdS2:IgdE2,IgdS:IgdE);
Vel.dtnum = TMEgridn(IgdS:IgdE);

ydaya    = datenum2yday(TMEgridn(IgdS));
Vel.yday = ydaya+Vel.dtnum-Vel.dtnum(1);

Vel.z = depthvec(IgdS2:IgdE2);
Vel.temp = ADP.temp(IgdS:IgdE);
Vel.xducer_depth = ADP.depth_xducer(IgdS:IgdE);
Vel.botdepth = InP.NomBotDepth;
Vel.lat = Lat;
Vel.lon = Lon;

Vel.MagDecUsed = MagDec;

Vel.info.snADCP = WHtyp(1:end);
Vel.info.MooringID = MooringID;
Vel.info.Cruise = CrName;
Vel.info.startyear = StartYear;
Vel.info.timeOffset = TmeOff;
Vel.info.timeOffset_notes = 'applied. Seconds that ADCP was ahead of GMT--linearly applied over timeseries';

Vel.info.sysconfig = sysconfig;
Vel.info.NomInstDepth = InP.NomInstDepth;
Vel.info.mfiles = sprintf(['in order:\n',...
                           'Conv_ADCP_mooring (coverts to mat files)\n'...
                           'Run_WHbeams_ALL.m (conv. to earth coord',...
                           ' & avgs pressure)\n']);
Vel.info.processing = ['10 1-s water pings & 10 BT pings avg. into ' num2str(TgridInt) ' minute ensembles.'];

if exist('Notes','var')==1;
  Vel.info.notes = Notes;
end




