% cruiseStruct_init_ADCP.m  -  initialize folders, indices for cruise;
% cd to the mfiles directory /cruisename/ADCP/mfiles before running this!
% Dave W, apr-2002, JBM Jul 2009, GV Mar 2015

% this is called by Conv_ADCP_mooring_CRUISENAME.m

% structures for matlab index files:
Cruise.name       = InP.Cruise.Name;
Cruise.start_date = InP.Cruise.StartTime;  % UTC, approx to start
Cruise.end_date   = InP.Cruise.EndTime;

[~,yrr]=datenum2yday(datenum(InP.Cruise.StartTime));

Def(1).start_year = yrr;
Def(1).start_dtnum = datenum(InP.Cruise.StartTime); 
Def(1).end_dtnum = datenum(InP.Cruise.EndTime); %or 
Def(1).dtnum_seconds_offset = 0;
Def(1).SWIMS_num = [];
Idx(1).dtnum_beg = [];
Idx(1).dtnum_end = [];
Idx(1).filename = [];
PROG ='cruiseStruct_init_ADCP.m';

%%%%%%%%%%%%%%%%%% HERE %%%%%%%%%%%%%%%%%%

data_typs={'ADCP'};

for id=1:length(data_typs)
  if id<=length(data_typs)
      fnidx = [data_typs{id} '_' InP.Cruise.Name '_matfiles.mat'];
    disp(['Init ' fnidx])
    if exist(fullfile('.',fnidx), 'file') && 0
      disp('   Already exists, skipping !!')
    else
      clear Set_params Index
      Set_params = Def;  Index = Idx;
      save( fullfile('.',fnidx),'Cruise','Set_params','Index','PROG');
    end
  end
end