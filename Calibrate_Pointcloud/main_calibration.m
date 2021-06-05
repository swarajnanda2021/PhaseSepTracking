%% Info
%-- main file for "Mechanical Calibration" using a calibration stage
%
%   $Author: Koen Muller 17-02-2018$
%
%-- modifications
%           name date what
%
%-- settings
%       [folder date (cal)]:  locate calibration files (Davis)
%       ctrl:               controls to processing
%       prop:               properties carried out experiments
%
%-- data management
%       .mat files  binary struct files: such as prop and ctrl
%       .dat files  binary data files for processing
%       .bin files  binary data size
%
%-- Todo
%
%	Check varslash implementation for linux
%
%-- notes

%% Matlab start
close all
clear all
clc

global folder date cal ctrl prop ext plotting Pmap Dmap

%% Additional Paths
dir= cd('..\'); % one folder up
addpath(genpath([cd '\Function_Libraries']))
cd(dir) %  back to original folder
clear dir

%% Locations
folder='D:\Desktop'; % Desktop Location
date='\2017_02_16'; % Date Folder
cal='\2017_02_16_Calibration'; % Calibration
test='\2017_02_16_Calibration'; % Test clicking

%% Controls
plotting='off'; % switch of all figures for cluster! ( 'off' | 'on')
[prop,ctrl] = proc_set; % processing setting function

%% def mappings
if 1
    [Pmap]=def_pmap; % projection map
    [Dmap]=def_dmap; % distortion map
end

%% load / process xml file
if 0
    proc_xmlfile; % here check that the reprojection is correct see comments
end

%% calibrate
if 0
    [Dpar,Kmat,Rmat,tvec,camprop]=calibrate;
end
% make camprop

%% camera warpings
if 0
    def_cmap
end

%% plot calibration
if 0
    plot_calpoints(Dpar,Kmat,Rmat,tvec);
end


%% perform simple example measurement
if 0
    click_raytracing(1) % input: frame
end

%% sandbox

