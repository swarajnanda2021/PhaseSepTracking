function [prop,ctrl,post,traj,eulr,grps,pair,mset] = proc_set
%proc_set Processing settings
% Output,
%   prop - Properties
%   ctrl - Control

%% global
global folder date rec

%% General
prop=importdata([folder date rec vsl 'Preproc3DTracking' vsl 'prop.mat']);
ctrl=importdata([folder date rec vsl 'Preproc3DTracking' vsl 'ctrl.mat']);
post=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'post.mat']);
traj=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'traj.mat']);
eulr=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'eulr.mat']);
grps=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'grps.mat']);
pair=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'pair.mat']);

%% movie controls
% time span movie
mset.tspan=[1 1000];%post.tproc;%[1 200];%

% ellipsoids
mset.ellipsoid = 'true';

% dynamic range and scaling
mset.cmod='minmax'; % nummeric | sqrt | log | imadjust | minmax
mset.drng=[200 2000]; % scaled | fixed by numeric: [100 1000]*prod(ctrl.rbox);

% time sequence
mset.tlen=10;
mset.lmin=15;
mset.tleng=1;

% domain and units
mset.dom=post.dom;
mset.col='normal'; % inverted [To be implemented]
mset.units='mm'; % object units

% movie and figure settings
mset.rend='-opengl'; % render setting, movie: none | figure: '-opengl' '-painters'
mset.ext='avi'; % extension, movie: 'avi' 'mp4' | figure: '-dbmp' '-depsc' 
mset.qual=100; % quality, movie: 0-100 | figure: '-r200' (... dpi)

%% adjust dynamic range
if isa(mset.drng,'double') % fixed nummeric
    switch mset.cmod
        case 'sqrt'
            mset.drng=sqrt(mset.drng);
        case 'log'
            mset.drng=log(mset.drng+1);
    end
    mset.drng=mset.drng*prod(ctrl.rbox);
end

%% save
save([folder date rec vsl 'Movie3DTracking' vsl 'mset.mat'],'mset'); % save controls
% save([folder date rec vsl 'Movie3DTracking' vsl 'fset'],'fset'); %

end