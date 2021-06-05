%% Info
%-- main file for processing movies
%
%   $Author: Koen Muller 14/06/19$
%
%-- settings
%       [folder date rec (cal)]:  locate carried out experiment (calibration)
%       ctrl:               controls to processing
%       prop:               properties carried out experiments
%       movie:              controls movie
%
%-- TODOs
%       Inverted colors image plane
%       Background removal image plane
%       add video extension option
%       add figure option? or maybe in conversion?

%% Matlab start
close all
clear
clc

global folder date rec cal plotting ctrl prop post traj eulr grps pair mset

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
%% create path
if exist([folder date rec vsl 'Movie3DTracking'],'dir')~=7
    mkdir([folder date rec vsl 'Movie3DTracking'])
end

%% Controls 
% dont show figure while processing
plotting='on'; 

% processing setting function
if 1
    [prop,ctrl,post,traj,eulr,grps,pair,mset]= proc_set;
else
    prop=importdata([folder date rec vsl 'Preproc3DTracking' vsl 'prop.mat']);
    ctrl=importdata([folder date rec vsl 'Preproc3DTracking' vsl 'ctrl.mat']);
    post=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'post.mat']);
    traj=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'traj.mat']);
    traj=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'locm.mat']);
    eulr=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'eulr.mat']);
    grps=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'grps.mat']);
    pair=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'pair.mat']);
    mset=importdata([folder date rec vsl 'Movie3DTracking' vsl 'mset.mat']);
end

%% Movie Camera Plane
if 0
    
    % Input by Name
    
    % Image data
    %   RawImage: Raw image plane
    %   DewarpedImage: Dewarped image plane
    %   BackgroundRemoval [under construction]
    
    % Camera data
    %   ImageEllipses: Ellipse identification
    %   ImageVelocimetry: Tracking in the image plane
    
    % Projected data
    %   ObjectEllipsoids ...
    %   ObjectVelocimetry: Project object tracking
    
    % Track singles
    %   ImageSingleEllipse: Follow one ellipse
    %   ObjectSingleEllipsoid: Project object tracking
    
%     !from plink  .. ImageSingleEllipse: Follow one ellipse
%     !separate windows ..  ObjectSingleEllipsoid: Project object tracking
%     !show index number
%     ! unify the above
    
    % Output by Name_Cam<N>
    movie_imageplane('RawImage')
    movie_imageplane('ImageEllipses')
    movie_imageplane('ObjectEllipsoids')
    movie_imageplane('ObjectVelocimetry')
	movie_imageplane('ObjectAndImageEllipsoids')
%    ! movie_imageplane('TrackingPerformance')
    
end

%% Movie Reference statistics
% if 0
%     movie_referencestats('TrackingPerformance') | ('LinearVelocity') | ('RotationRate')
% end

%% movie 3D tracking
if 0
    
    % Velocity
    movie_3dtracking('Velocity', 'Volume', 'TopView')
    movie_3dtracking('Velocity', 'Volume', 'Isometric')
    movie_3dtracking('Velocity', 'Volume', 'FrontView')
    movie_3dtracking('Velocity', 'Volume', 'CameraMotion')
    movie_3dtracking('Velocity', 'XYsliceZm2', 'CameraMotion')
    
    % Accelaration
    movie_3dtracking('Accelaration', 'Volume', 'TopView')
    
%    ! movie_3dtracking('TrackingPerformance', 'Volume', 'TopView')
    
end

%% Movie Eulerian Reference
if 0
    
    % Binned Velocity
    movie_eulerianreference('BinnedGraftieaux', 'GridNode', 'Volume', 'TopView')
    movie_eulerianreference('BinnedVelocity', 'GridNode', 'Volume', 'TopView')
    movie_eulerianreference('BinnedVelocity', 'GridNode', 'Volume', 'Isometric')
    movie_eulerianreference('BinnedVelocity', 'GridNode', 'Volume', 'FrontView')
    movie_eulerianreference('BinnedVelocity', 'GridNode', 'Volume', 'CameraMotion')
    movie_eulerianreference('BinnedVelocity', 'GridNode', 'XYsliceZ0', 'TopView')
    movie_eulerianreference('BinnedVelocity', 'GridPosition', 'XYsliceZm2', 'TopView')
    
    % Local polarization
    movie_eulerianreference('Polarization', 'GridNode', 'XYsliceZm2', 'TopView')
    
    % Linear Velocity / Magnitude
    movie_eulerianreference('LinearVelocity', 'GridNode', 'XYsliceZm2', 'TopView')  
    
    % Velocity
    movie_eulerianreference('Velocity', 'GridPosition', 'Volume', 'TopView')
    movie_eulerianreference('Velocity', 'GridPosition', 'Volume', 'FrontView')
    movie_eulerianreference('Velocity', 'GridPosition', 'Volume', 'CameraMotion')
    movie_eulerianreference('Velocity', 'GridPosition', 'XYsliceZm2', 'TopView')
    
    % Density
    movie_eulerianreference('Density', 'GridNode', 'XYsliceZm2', 'TopView')
    
    % Displacement (curved vector)
    % movie_eulerianreference('Displacement', 'GridPosition', 'Volume', 'TopView') % [slow & under construction]
    
    % FitResidual
    movie_eulerianreference('FitResidual', 'GridNode', 'XYsliceZm2', 'TopView') 
    
end

%% Movie spatial corrolation
if 0
    
    % Corrolation
    movie_spatialcorrelation('Correlation', 'XYsliceZ0', 'TopView')
    
    % Spectrum
    movie_spatialcorrelation('Spectrum', 'XYsliceZ0', 'TopView')
    
end

%% Movie Group motion
if 0
    
    % Input by name, 
    %   TopView_RotationAxisTrajectories 
    %   FrontView_RotationAxisTrajectories 
    %   TopView_FitVelocityField
    %   TopView_GroupVelocity
    %   TopView_GroupRotation
    %   TopView_GroupStrainRate
    
    % Output by name
    movie_groupmotion('RotationAxisTrajectories', 'Volume', 'TopView')
    movie_groupmotion('RotationAxisTrajectories', 'Volume', 'FrontView')
    movie_groupmotion('RotationAxisTrajectories', 'XYsliceZm2', 'TopView')
    
    % Velocity Field
    movie_groupmotion('FitVelocityField', 'Volume', 'TopView')
    movie_groupmotion('GroupVelocity', 'XYsliceZm2', 'TopView')
    movie_groupmotion('GroupRotation', 'XYsliceZm2', 'TopView')
    movie_groupmotion('GroupStrainRate', 'XYsliceZm2', 'TopView')
    
end

%% Movie Rotation Center

%% Movie Axis Rotation Center

%% Movie CotranslatingFrame

%% Movie CorotateFrame

%% Movie Lagrangian Reference[s?]
 % [under construction]

%% 3D trajectory frame
if 0
    % EllipsoidReconstruction | MidpointTracking
    movie_trajectframe('MidpointTracking') % [under construction]
    
    % Would be best to have reference measurement volume
    
end

%% 3D tangent frame
if 0
    % EllipsoidReconstruction | MidpointTracking
    movie_tangentframe('MidpointTracking') % [under construction]
    
    % Would be best to have reference measurement volume
    
end

%% Movie Visual Field
if 0
    % Projected Contours
    movie_visualfield('ProjectedContours') % [under construction]
    
    % Would be best to have reference measurement volume
    
end

%% Axisymetric coordinates
 % [under construction]

%% Sandbox

