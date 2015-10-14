function Vel = Conv_ADCP_mooring(ParamsFile)

% Vel = CONV_ADCP_MOORING(ParamsFile) - Convert set of ADCP files to raw
%   data files with lots of aux. data. Here ParamsFile is the params
%   setting m-file.. (leave off the '.m'!!).
%
%   Results going into a new set of matlab files.  Each RDI file will yield
%   a corresponding Matlab file.  The leading characters of the matlab
%   files are specified below, and are checked for prior use. The new index
%   file will also use this designation.  Access to these files is made via
%   get_ADCP_any.m. --  DPW, May 2008

% loading the parameters file (loads variable InP)
run(ParamsFile);

% For Alford mooring group, organize in ADCP s/n subfolders ...
mfiles = InP.mfiles; %#ok<*NODEF>

% Set root directory - everything is/will be under this.
% Cruise directory is a little bit misleading as it is really more a
% mooring directory in our case.
ADtop = fullfile([InP.Cruise.Dir,'/ADCP/']);
% eg /Users/gunnar/Projects/ttide/data/Moorings/A1/ADCP/

% Individual moored ADCPs, one for each deployment
ADnms = {['SN' InP.snADCP]};
DeployNo = InP.deployno; % for multi deploys of same S/N

ADsub = ADnms; % if multi deploys, each has own subfolder:
for i=1:length(DeployNo)
  if ~isnan(DeployNo(i))
    ADsub{i} = fullfile(ADnms{i},['Deployment' num2str(DeployNo(i))]);
  end
end
% eg 'SN22476_TTIDE'

Rnms = {[InP.snADCP]};
% eg '22476'

% Matlab root name for raw conversion files and index files
Mnms = {['SN' InP.snADCP '_' InP.Cruise.Name]};
% eg 'SN22476_TTIDE'

% Select deployment to convert to matlab in little gui:
[iA,ok] = listdlg('PromptString','Select an ADCP deployment:',...
    'Name','RDI-to-Matlab conversion', 'SelectionMode','single',...
    'ListString',ADsub);
if ~ok
  disp('None chosen, exit')
  return
end
% iA is ADCP/deployment index number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify files and folders %%
RawFld = fullfile(ADtop, ADsub{iA}, 'data_raw'); % raw BBsliced files
MatFld = fullfile(ADtop, ADsub{iA}, 'data_mat'); % matlab data files folder

RdiDeploy = Rnms{iA};

RawWild = ['*' InP.snADCP '*']; % to get BBsliced file names via 'dir'
MatRoot = Mnms{iA}; % for root names of matlab files and index
MatIndx = fullfile(MatFld, [MatRoot '_matfiles.mat']);

% load existing index, or initialize new one (THIS is for get_ADCP_any)
clear Cruise Set_params Index PROG
iSet = 1; % single-sequence index only
if exist(MatIndx,'file')
  load(MatIndx)
else
  % create the frigging index file! This actually creates Set_params
  cruiseStruct_init_ADCP
  Index(iSet).filename = [];
  Index(iSet).dtnum_beg = [];
  Index(iSet).dtnum_end = [];
  Index(iSet).RDIFile = [];
  Index(iSet).RDIbytes = [];
  PROG = mfiles;
end

% Gather RDI filenames into cell array (for non-empty files)
RDFs = dir(fullfile(RawFld,RawWild));
RDFnms = []; RDbyts = []; nf = 0;
for i=1:length(RDFs)
  if RDFs(i).bytes>10
    nf = nf+1;
    RDFnms{nf} = RDFs(i).name;
    RDbyts(nf) = RDFs(i).bytes;
  end
end
if nf<1
    disp( ['Exit, NO files found using: ' fullfile(RawFld,RawWild)] )
    return
end

RDcnv = 0*RDbyts; % conversion flag
% Search for new files to be converted, and for those that have grown
RDlist = [];
for i = 1:nf
  ia = find( strcmpi(RDFnms{i}, Index(iSet).RDIFile) );
  if isempty(ia)
    RDcnv(i) = 1; % new file
    RDlist{end+1} = [RDFnms{i} ' - new'];
  elseif RDbyts(i) > Index(iSet).RDIbytes(ia)
    RDcnv(i) = -ia; % File is larger than before
    RDlist{end+1} = [RDFnms{i} ' - GREW'];
  end
end

nf = length(find(RDcnv));
if nf<1
  disp(['Exit, no new/updated files found for ' fullfile(RawFld,RawWild)] )
  return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that processing is as intended %%
disp(['Looking in ' RawFld ','])
disp(['   found ' int2str(nf) ' RDI files:'])
disp('====================================')
disp(RDlist(:))
disp('====================================')
disp('The above-named files will be converted to matlab.')

% Create matlab data directory (may already exist)
if exist(MatFld,'dir')~=7
  [sta,msg] = mkdir( MatFld );
  if ~sta
    disp(msg);
    return
  end
end

disp('Resulting matlab files will be saved as (e.g.):')
disp( fullfile( MatFld, [MatRoot '*nn_00raw.mat'] ) )
disp([' with index = ' MatIndx])
disp(' ')
disp('Hit <enter> to continue:')
pause


for ifn = 1:length(RDcnv)
  if RDcnv(ifn)
    idxn = 0; % new file, add to end of index
    if RDcnv(ifn)<0
      idxn = -RDcnv(ifn); % existing, alter index entry
    end
    RDnm = RDFnms{ifn};
    typ = 2;  % THIS IS SAVING ALL DATA!! Set to <=0 to save less data

    %###########################%
    %### comment below if BT ###%
    %### or no.              ###%
    %###########################%
    
    %Vel = Get_ADCP_fullSC_BT( fullfile(RawFld,RDnm), -1, typ); %BT
    Vel = Get_ADCP_fullSC( fullfile(RawFld,RDnm), -1, typ); %NO BT

    % Apply time shift if necessary
    if isfield(InP,'TimeShift')==1; % time Shift is a constant shift
                                    % (e.g. in CRT..Canadian Retarded time,
                                    % NOT a drift over the timeseries).
      Vel.yday  = Vel.yday+InP.TimeShift./24;
      Vel.dtnum = Vel.dtnum+InP.TimeShift./24;
    end

    %checking read-in params vs. those in params file:
    if strcmp(Vel.coordXform,'Earth')==1;
      AA=1;
    elseif strcmp(Vel.coordXform,'Beam')==1;
      AA=0;
    elseif strcmp(Vel.coordXform,'Instrument')==1;
      AA=2;
    elseif strcmp(Vel.coordXform,'Ship')==1;
      AA=3;
    end

    % Sanity check coordinate xform
    if InP.set.Ptyp~=AA;
      disp(['warning, coordinate xform in params file does not match ',...
            'leader info in data file']);
      disp(['setting infile transform to ' Vel.coordXform ' coordinates']);
      disp('hit enter to continue');
      pause
    end
    % Set to correct value
    InP.set.Ptyp=AA;
    % Change the input parameter file to the correct xform settings.
    save(['SN' InP.snADCP '_' InP.MooringID '_params_' InP.Cruise.Name],'InP');

    % Check if there is any data in this file
    if isempty(Vel)
        disp('No data, skipping')
        continue
    end
    
    % Undo RDI depth calculation and use sw_dpth instead.
    % The RDI algorithm is the following:
    % Depth (dm) = Pressure(kPa) * (1.02-0.00069*ES)
    
    % ALSO FIXING WRAPPING PROBLEM: JM Aug 19, 2012
    % if any significant negative we have wrapping!!  ..or on the way up!
    % if exceed 3276.8 will wrap!! to negative values---starting at -3276.8
    AAA = Vel.depth_xducer;
    ibb = find(Vel.depth_xducer<0);
    if ~isempty(ibb); 
    AAA(ibb) = 3276.8.*2+AAA(ibb);
    %AA(ibb+1:iic(end))=AA(ibb+1:iic(end))+2.*abs(AA(ibb(1)+1));
    end
    ES    = Vel.params.ES; % estimated salinity (usually set to 35)
    prr   = AAA./(1.02-0.00069.*ES);
    DDnew = sw_dpth(prr',InP.Lat);
    Vel.depth_xducer = DDnew';
    
    % Beginning and end of time series
    ydb = Vel.dtnum(1); yde = Vel.dtnum(end);

    % Find subset that these data belong to (yearday range)
    % For a single mooring this is only one subset and it doesn't make
    % sense to check against the time input parameters and give and error
    % if we select the wrong time range. Why do we have to select this in
    % the input file anyways?
    iSet = 0;
    for i = 1:length(Set_params)
      if Set_params(i).start_dtnum <= ydb && Set_params(i).end_dtnum > ydb
        iSet = i;
        break
      end
    end
    if ~iSet
        disp(['Start dnum=' num2str(ydb) ...
            ' is outside the bounds of Grid Index - skipping!'])
        continue
    end

    % Create filename for raw .mat-file
    ia = length(RdiDeploy)+11;
    MTfn = [MatRoot RDnm(ia:end) 'raw.mat'];
    ib = strfind(MTfn, '.');
    if length(ib)>1
      MTfn(ib(1:end-1)) = '_';
    end
    
    % Append index numbers if more than one raw file
    if length(RDcnv)>1
      MTfn = [MTfn(1:end-4),'_',num2str(ifn),'.mat'];
    end
    
    % Save data in matlab file, update index
    disp(['  Saving as ' fullfile(MatFld, MTfn)])

    Vel.mfiles = mfiles;
    
    % Sometimes matlab will complain and not save your file if you don't
    % provide the -v7.3 switch...
    save(fullfile(MatFld, MTfn), 'Vel','-v7.3')

    % index
    if ~idxn
      idxn = length(Index(iSet).dtnum_end) + 1; % new entry
    end
    Index(iSet).filename{idxn}  = MTfn;
    Index(iSet).dtnum_beg(idxn) = ydb;
    Index(iSet).dtnum_end(idxn) = yde;
    Index(iSet).RDIFile{idxn}   = RDnm;
    Index(iSet).RDIbytes(idxn)  = RDbyts(ifn);
    
    
    save(MatIndx, 'Cruise','Set_params','Index','PROG')
  end
end

% Sort index by starting time
if iSet~=0;
  [~,iOr] = sort(Index(iSet).dtnum_beg);
  Index(iSet).filename  = Index(iSet).filename(iOr);
  Index(iSet).dtnum_beg = Index(iSet).dtnum_beg(iOr);
  Index(iSet).dtnum_end = Index(iSet).dtnum_end(iOr);
  Index(iSet).RDIFile   = Index(iSet).RDIFile(iOr);
  Index(iSet).RDIbytes  = Index(iSet).RDIbytes(iOr);
end

save(MatIndx, 'Cruise','Set_params','Index','PROG')

disp('All finished writing raw output.')


%reading data from all raw data files:

%once we've processed all of the raw files then we
%combine files if necessary after averaging,
%correct for declination
%grid onto a standard depth grid
%parse out bad data
%save into an informed directory structure.

%THEN combine ADCPs into one file if wanted.


