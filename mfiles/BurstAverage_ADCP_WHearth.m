function Vel = BurstAverage_ADCP_WHearth(VelAll,Filter)

% BurstAverage_ADCP_WHearth Average moored ADCP data collected in burst mode
%
%   Vel = BurstAverage_ADCP_WHearth(VelAll,Filter) will determine the burst
%   sampling pattern and average over all pings of one burst. The new time
%   vector is calculated from the time-mean over each burst.
%
%   INPUT   VelAll - Single ping ADCP data in earth coordinates (output
%                    from Run_WHbeams_ALL)
%           Filter (optional) - Set to >0 to filter out all pings that are
%                               outside Filter times one standard deviation
%                               of the burst mean.
%       
%   OUTPUT  Vel - Structure with velocities
%
%   Gunnar Voet
%   gvoet@ucsd.edu
%
%   Created: 2017-05-04

if nargin<2
  Filter = 0;
end

Vel = VelAll;

Time = Vel.dtnum;
DeltaTime = diff(Time);
t = Time(1:200);
dt = diff(t);

% Find the sampling frequency within the burst - this will be the dominant
% mode in the time diff
BurstPeriod = mode(DeltaTime);
fprintf(1,'time between pings in burst: %g seconds\n',BurstPeriod*24*3600);

% Now find the time between bursts. Need to subtract one BurstPeriod
ti = DeltaTime>2*BurstPeriod;
TimeBetweenBursts = mode(DeltaTime(ti))-BurstPeriod;
fprintf(1,'time between bursts: %g seconds\n',TimeBetweenBursts*24*3600);

% Find the first few times between bursts
tii = find(ti);
% Dominant mode of diff(tii) will give number of pings in one burst
BurstNumber = mode(diff(tii));

% % % last time stamp of burst
% % Time(tii(1))
% % % first time stamp of next burst
% % Time(tii(1)+1)

% start with the first full burst
tstart = Time(tii(1)+1);
tend = Time(tii(end));
tstarti = Vel.dtnum>=tstart & Vel.dtnum<=tend;
Vel = g_cutstruct(Vel,tstarti);


% Average matrix fields
MatrixFields = {'u','v','w','ecAve'};

for MI = 1:length(MatrixFields)

  % reshape for averaging
  ur = reshape(Vel.(MatrixFields{MI}),length(Vel.z),BurstNumber,[]);

  % averaging
  if ~Filter
    urm = squeeze(nansum(ur,2)./BurstNumber);
  else
    StandardDevation = nanstd(ur,1,2);
    StandardDevationMatrix = repmat(StandardDevation,1,10,1);
    MeanU = nanmean(ur,2);
    MeanUMatrix = repmat(MeanU,1,10,1);
    tmp = abs(ur-MeanUMatrix)-Filter.*StandardDevationMatrix;
    tmp2 = tmp;
    ur2 = ur;
    ur2(tmp2>0) = NaN;
    % tmp3 = reshape(tmp2,49,[]);
    urm = squeeze(nansum(ur2,2)./BurstNumber);
  end
  
  Vel.(MatrixFields{MI}) = urm;
  clear urm
  
end

% Average vectors
VectorFields = {'dtnum','yday','temp','xducer_depth'};

for MI = 1:length(VectorFields)

  % reshape for averaging
  ur = reshape(Vel.(VectorFields{MI}),BurstNumber,[]);

  % averaging
  urm = squeeze(nansum(ur,1)./BurstNumber);
  
  
  Vel.(VectorFields{MI}) = urm;
  clear urm
  
end