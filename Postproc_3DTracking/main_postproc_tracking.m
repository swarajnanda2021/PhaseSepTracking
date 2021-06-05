%% Info
%-- postprocessing file for "Lagrangian Particle Tracking"
%   
%   $Author: Koen Muller 29/03/2020$
%   
%-- data management
%       .mat files  binary struct files: such as prop and ctrl
%       .dat files  binary data files for processing
%       .bin files  binary data size
%       .h5 .csv files  final output for sharing
%   
%-- Notes
%   	These files maybe become a subfolder file tree that takes separate projects and plotting
%		Cleandata can improve on advanced segmentation e.g. frame wise segmentation criteria
%   
%-- Todo
%       All code have implemented variable slash for operating system but
%       need to be checked and where needed debugged

%% Matlab start
close all
clear all
clc

global folder date rec cal ctrl prop post traj locm eulr grps pair plotting

%% Additional Paths
dir= cd('..\'); % one folder up
addpath(genpath([cd '\Function_Libraries']))
cd(dir) %  back to original folder
clear dir

%% Locations
folder='C:\Work\PostGraduation\Scripts'; % Desktop Location
date='\ParticleTracking'; % Date Folder
rec='\Measurement37'; % Recording
cal='\Calibration_200514_141513'; % Calibration files
%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking'])
end

%% controls
plotting='off';
[prop,ctrl,post,traj,locm,eulr,grps,pair]  = proc_set;

%% initiate Tdata by cleaning Pobj
if 1
    postproc_cleandata % post : updated to choose longest tracks for duplicacies
end

%% kinematics from trajectories
if 1
    postproc_trajectories % traj
end

%% Track-quality analysis (also detects ghost track indices from inside the cavity)
if 1
   postproc_trackqualitymetric 
end


%% Spatial Mapping (binning into eulerian grid, removing all outliers)
if 1
    postproc_eulerianreference % eulr
end

%% Graftieaux vortex centred shifting (geometry based cavity velocity field averaging. centre located at grid resolution)
% also does cavity conditioned centering
if 1
    postproc_graftieaux 
end




% !%% tracking performance
% !if 0
%    ! postproc_trackingperformance % traj
% !end

%% frenet serret frame from trajectories
if 0
    postproc_frenetserret % [to be checked] traj
end

%% frequency statististics [waterfall]
if 0
    postproc_freqanalysis % [under construction] traj
end

% !%% autocorrelation
% !if 0
%    ! postproc_autocorrelation % traj
% !end

%% process reference statistics in local frame
if 0
    postproc_localmetrics % locm
end

%% process (revised version of) local metrics
% or maybe later to summarize what hass been analysed

%% define pairs in proximity
if 0
    postproc_pairs % pair , add mean velocity around object position
end

% !%% spatial structure
% !if 0
%    ! postproc_spatstructure % traj
% !end

%% separation analysis [finite time lyapunov]
if 0
    postproc_lyapanalysis % pair
end

%% Compute pair correlation
if 0
    postproc_paircorrelation % pair , add fluctuation velocity pair corrolation, add ranking to plot pair corrolation
end

%% process information transfer
if 0
%     postproc_!spatiotemporalcorrelation!then after that infotransfer %  [under construction] pair
end




%% Two Point Velocity Correlation
if 0
    postproc_spatialcorrelation %! its correlation % eulr
end

%% Process group frame
if 0
    postproc_groups % [add multiplegroups]  grps
end

%% Process instateneous center of rotation frame
if 0
    postproc_rotationcenterframe % grps
end

%% Process instateneoud axis of rotation frame
if 0
    postproc_rotationaxisframe % grps
end

%% Process co translating center of gravity frame
if 0
    postproc_cotranslatingframe % grps
end

%% Process co rotating center of gravity frame
if 0
    postproc_corotatingframe % [under construction, scale rotation to rigidly follow group] % grps
end

%% lagrangian frame to the group
if 0
    postproc_lagrangianreference % [under construction] % lagr
end

%% transform data to lagrangian frenetframe of the trajectory
if 0
    postproc_frenetframe % [under construction] % traj
end

%% transform data to lagrangian frame following the trajectory
if 0
    postproc_trajectframe % [under construction] traj
end

%% transform data to lagrangian frame of the trajectory
if 0
    postproc_tangentframe % [under construction] traj
end

%% transform trajectory frame to visual field
if 0
    postproc_visualfield % [error in azelaxis due missing toolbox] traj
end

%% process local axis symmetric coordinates
if 0
    postproc_axisymcoord % [error in azelaxis due missing toolbox] traj
end

%% process rank by different corrolation
if 0
    postproc_ranking % [under construction] locm
end

%% Sandbox
% linear system identifycation
