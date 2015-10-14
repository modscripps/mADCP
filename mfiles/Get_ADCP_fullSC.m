function Vel = Get_ADCP_fullSC(RDIname, filetype, savetype, ens_ind_rng)

% Vel = Get_ADCP_fullSC(RDIname, filetype, savetype, ens_ind_rng)
%
%	Read in RDI file(s), put requested results in structure = Vel
%
%  INPUT
%  filetype: 0=process single file=RDIname, 1=multi LTA files
%		starting with RDIname (more later);
%	 savetype: 0=just u,v,w, depth,time,position vectors,
%		1=more stuff (%good, etc), 2=everything
%  ens_ind_rng: min,max index values for ensembles in raw RDI file to
%       parse and convert to matlab data (11/04)
% 
%  OUTPUT
%  Vel = structure with useful RDI data
%
%  Dave W, 03-Apr-2001, 16-nov-2004
%  Get_ADCP_profRaw.m modified for raw ping data received from SWIMS ADCPs
%
%  GV making a few changes here and there, mostly cleaning up the code



% GV: comment out the following, not needed here
% global InP

if nargin<2
   filetype = 0; % single profile file (EG via serial Ensemble output)
end
if nargin<3
   savetype = 0;
end
if nargin<4
    ens_ind_rng = []; % get all ensembles
end

% filetype = -1; % for Raw: single file with multiple ensembles. A moored
% ADCP would be filetype = -1

Vel = [];

if filetype==0  % check that this is one complete file/profile
   fid=fopen(RDIname, 'r','ieee-le');
   if fid<0
      error(['Cannot open file = ' RDIname])
   end
   A = fread(fid,'char');
   frewind(fid);
   id = fread(fid,1,'uchar');
   src= fread(fid,1,'uchar');
   nbytes=fread(fid,1,'ushort');
   fclose(fid);
   if id~=127 || src~=127
      error([RDIname ' does not start with RDI header record']);
   end
   if length(A)-nbytes ~= 2
      error([RDIname ' is wrong size for single RDI profile']);
   end
end % of solo-file check

% raw2mat save to a file. Create an arbitrary filename based on current
% timestamp.
[p,n,e] = fileparts(RDIname);
x = datevec(now); x = round(x(6));
TMPfil = ['TMP' num2str(x) '_' n];

% GV: why are these global? this is strange and shouldn't be needed
global year month day hour minute second

% Convert File(s)
if filetype < 1 % either: 0 = solo profile, or <0 = one file, multi-profs
  raw2mat_SCv3(RDIname,TMPfil, [], {'.mat'}, [], ens_ind_rng);
else
  error('Get_ADCP_fullSC.m only works for filetype<0')
end

year=[];month=[];day=[];hour=[];minute=[];second=[];

% load in Matlab arrays, convert to useful units, save in structure
load(TMPfil)
delete([TMPfil '.mat']);

% Check if no valid data were found, exit if so
if isempty(ensemble_number)
    Vel = [];
    return
end

BADVEL = -32768;


%interpreting sysconfig

Oxx=dec2bin(sysconfig); %refer to table in ADCP output format

LSB=Oxx(end-7:end);
MSB=Oxx(1:end-8);

if str2num(LSB(1))==1;
  beamdir='upfacing';
else
  beamdir='downfacing';
end

if str2num(LSB(5))==1;
  beampattern='convex';
else
  beampattern='concave';
end

if str2num(LSB(6:end))==0;
  adcpFreq='75kHz';
elseif str2num(LSB(6:end))==1;
  adcpFreq='150kHz';
elseif str2num(LSB(6:end))==10;
  adcpFreq='300kHz';
elseif str2num(LSB(6:end))==11;
  adcpFreq='600kHz';   
elseif str2num(LSB(6:end))==100;
  adcpFreq='1200kHz'; 
elseif str2num(LSB(6:end))==101;
  adcpFreq='2400kHz'; 
end

% MSB may not be 8 characters long, need to pad
% with zeros if this is not the case.
if str2num(MSB(end-1:end))==0;
  beamAngle='15';
elseif str2num(MSB(end-1:end))==1;
  beamAngle='20';
elseif str2num(MSB(end-1:end))==10;
  beamAngle='30';
elseif str2num(MSB(end-1:end))==11;
  beamAngle='nonstandard';
end

UpDown = bitand(sysconfig,2^7);
UpDown = dec2bin(UpDown);
UpDown = UpDown(1);

if UpDown % up-facing bit == 1
  UpDown = 1; % -1=down, +1=upward
else
  UpDown = -1;
end

%interpreting beam coordinates
Cxx = dec2bin(coords); %refer to table in ADCP output format

%LSB=Cxx(end-7:end);
%MSB=Oxx(1:end-8);
if length(Cxx)==1;
  XA=str2num(Cxx);
  if XA==0;
    coordXform='Beam';
  elseif XA==1;
    coordXform='Instrument';
  else
    disp('confused coordinate system. Crashing.');
    return
  end
elseif length(Cxx)>=5;

if str2num(Cxx(4:5))==11;
  coordXform='Earth';
elseif str2num(Cxx(4:5))==1;
  coordXform='Instrument';
elseif str2num(Cxx(4:5))==10;
  coordXform='Ship';
elseif str2num(Cxx(4:5))==0;
  coordXform='Beam';
end

else
  disp('confused coordinate system. Crashing.');
  return
end

Vel.ens_no = ensemble_number + ensMSB*65536;

% Some will be the same for all ensembles:
Vel.z_adcp =   UpDown*(dis1/100 + [0:nbins-1]*binlen/100); % meters  %this should be negative up, positive down.


%Vel.p_adcp = Vel.z_adcp; % db

%Vel.headbias=headbias;




Vel.sysconfig.beamdir     = beamdir;
Vel.sysconfig.beampattern = beampattern;
Vel.sysconfig.freq        = adcpFreq;
Vel.sysconfig.beamangle   = beamAngle;

%setting params;
Vel.params.WX = pulselen./100;  %transmit pulse length / m
Vel.params.WC = lowcorr;  %low correlation threshold
Vel.params.WE = errthresh;   % velocity error threshold
Vel.params.WG = goodthresh;   %min. percent good
Vel.params.EB = headbias;  %heading bias
Vel.params.ES = estsal(1);


if waterband==0;
  waterband = 'broad';
else
  waterband = 'narrow';
end
Vel.params.WB = waterband;       %bandwidth: broadband=0, narrow=1;

Vel.params.names={'WX=pulselength /m';'WC=low corr. threshold';'WE=error threshold';'WG=perct. good threshold';...
    'EB=heading bias';'ES=estimated salinity';'WB=bandwidth mode'};

% Others will vary:
year=year+2000;

%Vel.yday = yearday(day,month,year,hour,minute,second+hundreths/100);
ydayS=yearday(day(1),month(1),year(1),hour(1),minute(1),second(1)+hundreths(1)/100);
ydayEL=datenum(year,month,day,hour,minute,second+hundreths./100)-datenum(year(1),month(1),day(1),hour(1),minute(1),second(1)+hundreths(1)./100);

Vel.yday=ydayS+ydayEL;
Vel.dtnum=datenum(year,month,day,hour,minute,second+hundreths./100);
% %% BT - range and velocity
% if isempty(btrange1) % If set for NOT bottom-tracking, put in NaNs (3-18-2003)
%         ix = NaN*Vel.ens_no;
%         btrange1 = ix; btrange2 = ix; btrange3 = ix; btrange4 = ix;
%         btvel1 = ix; btvel2 = ix; btvel3 = ix; btvel4 = ix;
% end

% Orientation data (Internal for SWIMS)
Vel.heading = heading/100;
Vel.pitch = pitch/100;
Vel.roll = roll/100;
Vel.hdg_std = stdhed/1;
Vel.pitch_std = stdpitch/10;
Vel.roll_std = stdroll/10;
Vel.depth_xducer = xducerdepth/10;

%FIXING WRAPPING PROBLEM: JM Aug 19, 2012
%ibb=find(diff(Vel.depth_xducer)<-6000);
%if any significant negative we have wrapping!!  ..or on the way up!

%if exceed 3276.8 will wrap!! to negative values---starting at -3276.8

%ibb=find(diff(Vel.depth_xducer)<-100);
AA=Vel.depth_xducer;
ibb=find(Vel.depth_xducer<0);
if ~isempty(ibb); 

AA(ibb)=3276.8.*2+AA(ibb);
%AA(ibb+1:iic(end))=AA(ibb+1:iic(end))+2.*abs(AA(ibb(1)+1));
%Vel.depth_xducer=AA;
end

% GV: moved this into Conv_ADCP_mooring_[experiment].m to have latitude
% available. Not good to hardcode latitude here!

% correcting depth based on better algorithm. The RDI algorithm is the
% following:
% Depth (dm) = Pressure(kPa) * (1.02-0.00069*ES),  
%correcting depth calculation using seawater toolbox.
% % % ES  = Vel.params.ES; %estimated salinity.
% % % prr = AA./(1.02-0.00069.*ES); %should get ES (usually 35)
% % % %DDnew=sw_dpth(prr',InP.Lat);
% % % DDnew=sw_dpth(prr',71);
% % % Vel.depth_xducer=DDnew';

Vel.soundvel   = soundspeedRDI;
Vel.temp       = degC/100;
Vel.coordXform = coordXform;


% add A/D channels - MHA 8/12/2012
if exist('adcchan0','var')==1
Vel.adcchan0 = adcchan0;
end
if exist('adcchan1','var')==1
Vel.adcchan1 = adcchan1;
end



%% Measured beam components
bad = find(v1==BADVEL); if bad, v1(bad) = NaN; end
bad = find(v2==BADVEL); if bad, v2(bad) = NaN; end
bad = find(v3==BADVEL); if bad, v3(bad) = NaN; end
bad = find(v4==BADVEL); if bad, v4(bad) = NaN; end

% These are in earth's frame for typical SC deployment, but may be beam vels:
if strcmp(coordXform,'Beam')==1;
Vel.v1 = v1/1000; Vel.v2 = v2/1000; % m/s
Vel.v3 = v3/1000; Vel.v4 = v4/1000; % m/s
else
Vel.u = v1/1000; Vel.v = v2/1000; % m/s
Vel.w = v3/1000; Vel.werr = v4/1000; % m/s
end


        
Vel.v_info={'v1,v2,v3,v4 are beams 1-4 in beams coords,';...
    'and u,v,w,werr in earth coords.'};
if savetype>0 
    % Save more ...
    Vel.ec1_bm = e1; Vel.ec2_bm = e2; Vel.ec3_bm = e3; Vel.ec4_bm = e4;
    Vel.ec_info={'echo intensity at about 0.45 db per count for each beam'};
    
    
    
    Vel.cor1_bm = cor1; Vel.cor2_bm = cor2; Vel.cor3_bm = cor3; Vel.cor4_bm = cor4;
    Vel.cor_info={'correlation magnitude for each beam, 0-255, 255=perfect corr.'};
    
    if ~isempty(find(~isnan([pg1(:); pg2(:); pg3(:); pg4(:)])));  %this will be empty when conversion from beam coords.
        Vel.pg1 = pg1; Vel.pg2 = pg2; Vel.pg3 = pg3; Vel.pg4 = pg4;
        Vel.pg_info={'% of good data collected';'when in beam coords. % of each beam';...
            '     when in earth coords';...
            'pg1=% of 3-beam transformations because WC too low';'pg2=% of transformations rejected because WE exceeded';...
            'pg3=% of more than 1 beam bad in bin';'pg4=% of 4-beam transformations'}; 
    end
end


Vel.more_info={'z_adcp: negative upwards';...
    'roll,pitch,heading in deg'};


return
