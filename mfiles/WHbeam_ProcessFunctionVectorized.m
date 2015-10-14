function ADP = WHbeam_ProcessFunctionVectorized(ADRaw,UpDown,ZGrid,theta_o,Cnvx,MagDec,WC_val,ExcludeBins)

% WHbeam_ProcessFunctionVectorized (well, partly, still looping a lot)
%
% This used to be the script with the silly global variable UpDown

% Check for exclude bins
if nargin<8
  ExcludeBins = 0;
end

% initialize output structure
ADP.depth = ZGrid; % surface-relative depth grid
Trow = length(ADP.depth);
Tcol = length(ADRaw.dtnum);

ADP.u = NaN * ones(Trow,Tcol);
ADP.v = ADP.u;
ADP.w = ADP.u;
ADP.werr = ADP.u;

ADP.ec1 = 0*ones(Trow,Tcol);
ADP.ec2 = ADP.ec1;
ADP.ec3 = ADP.ec1;
ADP.ec4 = ADP.ec1;

ADP.dtnum = ADRaw.dtnum;
ADP.ens_no = ADRaw.ens_no;
ADP.soundvel = ADRaw.soundvel;
ADP.svel_calc = ADP.soundvel;

if UpDown==0 %up facing
ADP.z_adcp = -ADRaw.z_adcp;
else
ADP.z_adcp = ADRaw.z_adcp; 
end

ADP.depth_xducer = ADRaw.depth_xducer;
vvs = {'svel_calc','ens_up','yd_up','temp'};
for i=1:length(vvs)
    if isfield(ADRaw, vvs{i}) && ~isempty(ADRaw.(vvs{i}))
        ADP.(vvs{i}) = ADRaw.(vvs{i});
    end
end

if size(ADP.z_adcp,1)==1
    ADP.z_adcp = ADP.z_adcp'; % column vector
end

% WC_val=10; % GV: used to be hardcoded here, now as argument for the
% function

% Screen data for 'gross' problems
ix = find(ADRaw.cor1_bm(:) < WC_val);
ADRaw.v1(ix) = NaN;
ix = find(ADRaw.cor2_bm(:) < WC_val);
ADRaw.v2(ix) = NaN;
ix = find(ADRaw.cor3_bm(:) < WC_val);
ADRaw.v3(ix) = NaN;
ix = find(ADRaw.cor4_bm(:) < WC_val);
ADRaw.v4(ix) = NaN;

% RDI transformation matrix from beams 1-4 to
%   u(1-2), v(4-3), w(avg xz,yz), err vel 
Bm2InTx = beam2inst(theta_o, Cnvx);
% Make range vector from the depth vector (before pitch/roll adjustments)
c_tho = cos(theta_o*pi/180);
s_tho = sin(theta_o*pi/180);
rBM_o = ADP.z_adcp/c_tho;
% soundspeed correction, actual/nominal
if isempty(ADP.svel_calc)
    ADP.svel_calc = ADP.soundvel;
end
SSadj = ADP.svel_calc ./ ADP.soundvel; 
% % Turn the recorded pitch into real pitch... insignif. for small rolls.
% ADRaw.pitch=180/pi*atan(tan(pi/180*ADRaw.pitch).*cos(pi/180*ADRaw.roll));
% % Slighty filter pitch and roll:
% ADP.pitch = [ ADRaw.pitch(1), ...
%         (0.2*ADRaw.pitch(1:end-2)) + (0.6*ADRaw.pitch(2:end-1)) + (0.2*ADRaw.pitch(3:end)), ...
%         ADRaw.pitch(end) ];
% ADP.roll = [ ADRaw.roll(1), ...
%         (0.2*ADRaw.roll(1:end-2)) + (0.6*ADRaw.roll(2:end-1)) + (0.2*ADRaw.roll(3:end)), ...
%         ADRaw.roll(end) ];
%% Just copy pitch and roll (and compass heading) instead:
ADP.pitch = ADRaw.pitch;
ADP.roll = ADRaw.roll;
ADP.heading = ADRaw.heading;
% Signs for correcting range(depth)
if UpDown == 1
    sg1 = 1; sg3 = 1;
else
    sg1 = 1; sg3 = -1;
end


%% Vectorize the transformation to avoid looping

% The first step is to compute depth vectors for each beam.  These differ from 
% each other since the beams are tilted; Then adjust for sound speed.
% Roll affects beams 1 and 2,  pitch affects beams 3 and 4.
rBM = rBM_o; % adjust depths for pitched/rolled beam angles
rBMV = repmat(rBM(:),1,length(ADP.dtnum));
rollV = repmat(ADP.roll(:)',length(rBM),1);
pitchV = repmat(ADP.pitch(:)',length(rBM),1);
z1_rel = cos( pi/180.*(theta_o +  sg1.*rollV ) ).*rBMV;
z2_rel = cos( pi/180.*(theta_o + -sg1.*rollV ) ).*rBMV; 
z3_rel = cos( pi/180.*(theta_o +  sg3.*pitchV) ).*rBMV;
z4_rel = cos( pi/180.*(theta_o + -sg3.*pitchV) ).*rBMV;

% compute absolute depths considering sound speed and ADCP depth
z1_abs = z1_rel;
z2_abs = z2_rel;
z3_abs = z3_rel;
z4_abs = z4_rel;
z1_abs(1,:) = z1_rel(1,:) .* SSadj + ADP.depth_xducer;
z1_abs(2:end,:) = repmat(z1_abs(1,:),length(rBM)-1,1) + cumsum(diff(z1_rel,1,1).*repmat(SSadj,length(rBM)-1,1),1);
z2_abs(1,:) = z2_rel(1,:) .* SSadj + ADP.depth_xducer;
z2_abs(2:end,:) = repmat(z2_abs(1,:),length(rBM)-1,1) + cumsum(diff(z2_rel,1,1).*repmat(SSadj,length(rBM)-1,1),1);
z3_abs(1,:) = z3_rel(1,:) .* SSadj + ADP.depth_xducer;
z3_abs(2:end,:) = repmat(z3_abs(1,:),length(rBM)-1,1) + cumsum(diff(z3_rel,1,1).*repmat(SSadj,length(rBM)-1,1),1);
z4_abs(1,:) = z4_rel(1,:) .* SSadj + ADP.depth_xducer;
z4_abs(2:end,:) = repmat(z4_abs(1,:),length(rBM)-1,1) + cumsum(diff(z4_rel,1,1).*repmat(SSadj,length(rBM)-1,1),1);

% map beam velocities onto standard depths (sfc-refn'd),
%   adjust for actual soundspeed
u1 = nan(length(ADP.depth),length(ADP.dtnum));
u2 = u1;
u3 = u1;
u4 = u1;

% for ci = 1:length(ADP.dtnum)
%   u1(:,ci) = interp1(z1_abs(:,ci),ADRaw.v1(:,ci),ADP.depth).*SSadj(ci);
%   u2(:,ci) = interp1(z2_abs(:,ci),ADRaw.v2(:,ci),ADP.depth).*SSadj(ci);
%   u3(:,ci) = interp1(z3_abs(:,ci),ADRaw.v3(:,ci),ADP.depth).*SSadj(ci);
%   u4(:,ci) = interp1(z4_abs(:,ci),ADRaw.v4(:,ci),ADP.depth).*SSadj(ci);
% end

% Look for bins to exclude
zi = 1:length(ADRaw.z_adcp);
zi = logical(ones(1,length(ADRaw.z_adcp)));
if ExcludeBins
  fprintf(1,'WHbeam_ProcessFunctionVectorized: excluding bin %02d in interpolation\n',ExcludeBins);
  zi(ExcludeBins) = 0;
  fprintf(1,'\n');
end

% use faster interp1 (but need to have monotonically increasing x):
ztmp = nanmean(z1_abs,2);
if ztmp(1)>ztmp(end)
z1_abs_flud = flipud(z1_abs);
v1flud = flipud(ADRaw.v1);
z2_abs_flud = flipud(z2_abs);
v2flud = flipud(ADRaw.v2);
z3_abs_flud = flipud(z3_abs);
v3flud = flipud(ADRaw.v3);
z4_abs_flud = flipud(z4_abs);
v4flud = flipud(ADRaw.v4);
zi = fliplr(zi);
else
z1_abs_flud = z1_abs;
v1flud = ADRaw.v1;
z2_abs_flud = z2_abs;
v2flud = ADRaw.v2;
z3_abs_flud = z3_abs;
v3flud = ADRaw.v3;
z4_abs_flud = z4_abs;
v4flud = ADRaw.v4;
end

% Progress bar
fprintf('beam2earth: interpolating to absolute depth grid\n');
dtn = length(ADP.dtnum);
dtprct = round(prctile(1:dtn,10:10:90));
fprintf('|->%s|  0%%',repmat(' ',1,21))

% Interpolate
for ci = 1:length(ADP.dtnum)
  [lia,locb] = ismember(ci,dtprct);
  if lia
    fprintf(repmat('\b',1,28))
    fprintf('|%s->%s| %d%%',repmat('-',1,locb*2),repmat(' ',(10-locb)*2,1),locb.*10)
  end
  u1(:,ci) = interp1qr(z1_abs_flud(zi,ci),v1flud(zi,ci),ADP.depth).*SSadj(ci);
  u2(:,ci) = interp1qr(z2_abs_flud(zi,ci),v2flud(zi,ci),ADP.depth).*SSadj(ci);
  u3(:,ci) = interp1qr(z3_abs_flud(zi,ci),v3flud(zi,ci),ADP.depth).*SSadj(ci);
  u4(:,ci) = interp1qr(z4_abs_flud(zi,ci),v4flud(zi,ci),ADP.depth).*SSadj(ci);
end
fprintf(repmat('\b',1,28))
fprintf('|%s->|100%%\n',repmat('-',1,21))

% GV: I couldn't find a way to vectorize this - maybe there is none as
% every column in the vector has different y values...


% Tried the following, but this also interpolates along x which we do not
% want...
% Y = z1_abs;
% Yr = reshape(Y,[],1);
% X = repmat(ADP.dtnum,length(rBM),1);
% Xr = reshape(X,[],1);
% Z = ADRaw.v1;
% Zr = reshape(Z,[],1);
% Yq = repmat(ADP.depth(:),1,length(ADP.dtnum));
% Yqr = reshape(Yq,[],1);
% Xq = repmat(ADP.dtnum,length(ADP.depth),1);
% Xqr = reshape(Xq,[],1);
% 
% 
% F = scatteredInterpolant(Xr,Yr,Zr,'linear','none');
% u1f = F(Xqr,Yqr);
% u1fr = reshape(u1f,length(ADP.depth),[]);

ADP.ec1 = nan(length(ADP.depth),length(ADRaw.dtnum));
ADP.ec2 = ADP.ec1;
ADP.ec3 = ADP.ec1;
ADP.ec4 = ADP.ec1;

if ztmp(1)>ztmp(end)
ec1flud = flipud(ADRaw.ec1_bm);
ec2flud = flipud(ADRaw.ec2_bm);
ec3flud = flipud(ADRaw.ec3_bm);
ec4flud = flipud(ADRaw.ec4_bm);
else
ec1flud = ADRaw.ec1_bm;
ec2flud = ADRaw.ec2_bm;
ec3flud = ADRaw.ec3_bm;
ec4flud = ADRaw.ec4_bm;
end


% map echo intensities, for later screening
fprintf('beam2earth: map echo intensities\n');
fprintf('|->%s|  0%%',repmat(' ',1,21))
for ic = 1:length(ADP.dtnum)
  
  [lia,locb] = ismember(ic,dtprct);
  if lia
    fprintf(repmat('\b',1,28))
    fprintf('|%s->%s| %d%%',repmat('-',1,locb*2),repmat(' ',(10-locb)*2,1),locb.*10)
  end
  
ADP.ec1(:,ic) = interp1qr( z1_abs_flud(:,ic), ec1flud(:,ic), ADP.depth);
ADP.ec2(:,ic) = interp1qr( z2_abs_flud(:,ic), ec2flud(:,ic), ADP.depth);
ADP.ec3(:,ic) = interp1qr( z3_abs_flud(:,ic), ec3flud(:,ic), ADP.depth);
ADP.ec4(:,ic) = interp1qr( z4_abs_flud(:,ic), ec4flud(:,ic), ADP.depth);
end
fprintf(repmat('\b',1,28))
fprintf('|%s->|100%%\n',repmat('-',1,21))

ADP.ec1 = round(ADP.ec1);
ADP.ec2 = round(ADP.ec2);
ADP.ec3 = round(ADP.ec3);
ADP.ec4 = round(ADP.ec4);


% Now, transform beam velocities to earth (geomagnetic) coordinates.
fprintf('beam2earth: transform beam to earth coordinates\n');
fprintf('|->%s|  0%%',repmat(' ',1,21))
for ic = 1:length(ADP.dtnum)
  
  [lia,locb] = ismember(ic,dtprct);
  if lia
%     fprintf('%d%% ',locb.*10)
    fprintf(repmat('\b',1,28))
    fprintf('|%s->%s| %d%%',repmat('-',1,locb*2),repmat(' ',(10-locb)*2,1),locb.*10)
  end
  
In2Geo = inst2earth(ADP.heading(ic)+MagDec, ADP.pitch(ic), ADP.roll(ic), UpDown);
VelsGeo = ( In2Geo * Bm2InTx * [u1(:,ic) u2(:,ic) u3(:,ic) u4(:,ic)]' )';
    ADP.u(:,ic) = VelsGeo(:,1);
    ADP.v(:,ic) = VelsGeo(:,2);
    ADP.w(:,ic) = VelsGeo(:,3);
    ADP.werr(:,ic) = VelsGeo(:,4);
end
fprintf(repmat('\b',1,28))
fprintf('|%s->|100%%\n',repmat('-',1,21))
fprintf('\n')
%     VelsGeo = ( In2Geo * Bm2InTx * [u1 u2 u3 u4]' )';
%     % Gather components into output structure
%     ADP.u(:,ic) = VelsGeo(:,1);
%     ADP.v(:,ic) = VelsGeo(:,2);
%     ADP.w(:,ic) = VelsGeo(:,3);
%     ADP.werr(:,ic) = VelsGeo(:,4);


return
