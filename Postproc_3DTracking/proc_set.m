function [prop,ctrl,post,traj,locm,eulr,grps,pair] = proc_set
%proc_set Processing settings
% Output,
%   prop - Properties
%   ctrl - Control
%   +++

%% global
global folder date rec 

%% Load
prop=importdata([folder date rec vsl 'Preproc3DTracking' vsl 'prop.mat']);
ctrl=importdata([folder date rec vsl 'Preproc3DTracking' vsl 'ctrl.mat']);

%% Controls



%-- #) post processing and physical segmentation
post.tproc=[1 1000];%[1 100];%
post.dom=[-30 30
    -30 30
    -10 40]; % domain size [x [min max] ; y ; z  ]
post.tlen=1; % track length
post.cam=inf; % segment minimum camera cover tracks (0-all to inf-none)
post.rpe=1; % segment reprojection error (0-all to inf-none)
post.vel=inf; % segment velocity (0-all to inf-none)
post.siz=inf; % segment size quadrics (0-all to inf-none)
prop.modelname = 'Training_08122020_batchsize96_epochs50_ResNet18_titanx_optionset_2'; % Make sure you put the same model you used for inference

% post.grid=[-2 0.2 2
%     -2 0.2 2
%     0.3 0.2 2.3];

%-- #) trajectory shape
traj.tfit=5; % number time-frames 
traj.tord=3; % time poly-order
traj.outl='none'; % dynamic outlier filter ( none | median | meanstd | quantile | prctile | bimodal )
traj.conf=2; % outlier confidence
traj.iter=3; % number iterations (1 - inf)
%   traj.freq

%-- #) local metrics
locm.nngb=20; % number of [ranking] neighbors

%-- #) eulerian reference
eulr.gres=2; % grid resolution in units obj space
eulr.sres=2; % 0.5; % spatial resolution [coherent motion] in units obj space
eulr.tres=5; % temporal resolution [coherent motion] window in frms
eulr.outl='meanstd'; % outlier filter ( none | median | meanstd | quantile | prctile | bimodal )
eulr.conf=3; % outlier confidence (here by uniform distribution) for meanstd
% eulr.conf = [25 75];
eulr.iter=3; % number iterations (1 - inf)
eulr.lmin=5;
% Lagrangian and frame instateneous
%   lagr

%-- #) groupframe
grps.methods='nmodal'; % hartigan test
traj.tord=3; % time poly-order
grps.tres=7; % temporal resolution [coherent motion] window in frms
grps.outl='median'; % outlier filter ( none | median | meanstd | quantile | prctile | bimodal )
grps.conf=5; % outlier confidence (here by uniform distribution)
grps.iter=3; % number iterations (1 - inf)

%-- #) information transfer
pair.mran=20; % range in m
pair.cwin=10; % corrolation window size and intersection
pair.pdis=1; % maximum pair promiximity distance (to define s possible pair)
pair.pfrm=15; % minimum pair promiximity time frames (to define a possible pair)
% pair.sampl=25; % pair sampling for binning the data over unit domains
% pair.lyap

%% save proc set
disp('save properties and controls')
save([folder date rec vsl 'Postproc3DTracking' vsl 'post'],'post'); % save controls
save([folder date rec vsl 'Postproc3DTracking' vsl 'traj'],'traj'); % save controls
save([folder date rec vsl 'Postproc3DTracking' vsl 'locm'],'locm'); % save controls
save([folder date rec vsl 'Postproc3DTracking' vsl 'eulr'],'eulr'); % save controls
save([folder date rec vsl 'Postproc3DTracking' vsl 'grps'],'grps'); % save controls
save([folder date rec vsl 'Postproc3DTracking' vsl 'pair'],'pair'); % save controls
% save([folder date rec vsl 'Postproc' vsl 'lagr'],'post'); % save controls

end