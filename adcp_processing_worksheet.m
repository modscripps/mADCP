% Run processing steps for a sample moored ADCP
%
% Gunnar Voet
% gvoet@ucsd.edu
%
% Created: 09/29/2015

clear

% Add path for adcp processing toolbox
addpath ADCP_Processing/
BaseDirectory = 'sample_data/';

Project = 'TEST';


%% AM1 WH 3160 - DONE
sn  = '3160';
mid = 'AM1';

edit set_3160_paramsTEST.m

Conv_ADCP_mooring(sprintf('set_%s_params%s',sn,Project))

% All data
VelAll = Run_WHbeams_ALL([],[],sprintf('SN%s',sn),mid,Project);
save(fullfile(BaseDirectory,...
     sprintf('ADCP/SN%s/data_mat/SN%s_%s_ALL.mat',sn,sn,Project)),'VelAll')

% Time average 10min
Vel = Average_ADCP_WHearth(VelAll, 1, 1);
save(fullfile(BaseDirectory,...
  sprintf('ADCP/SN%s/data_mat/SN%s_%s_AVE_10min.mat',sn,sn,Project)),'Vel')