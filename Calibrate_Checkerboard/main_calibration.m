%% Info
%   
%   Author: Koen Muller, 2019
%   
%   N-view camera calibration in separate consecutive parts:
%       Part I: (reproducing include files)
%           Image processing
%
%       Part II: (reproduce include files)
%           Distortion correction rectifying the images
%           Single view calibration
%           Consistent multiple view calibration
%           Nummeric evaluation of camera mapping
%
%       Part III: (perform physical measurement)
%           Several functions to evaluate the results that can be modified
%       	internally
%   
%   The if 0 / 1 statement can be modified to run a certain part or
%   configuration of scripts by the user.
%
%   Data formats
%       .mat files matlab files
%       .dat files unstructured binary data files
%       .bin files corresponding binary data size to .dat file
%
%   Notes
%       i) Camera coordinate system is define as [x y k] where x and y are
%           related to the image plane and k is the depth of field. This
%           coordinate system is therefore left-handed.
%       ii) The world coordinates is defined by [X Y Z] where X is the 
%           width, Y is the (visual) depth and Z is the height of the 
%           measurement volume. This coordinate system is right-handed.

%!restructure output files tree
%!export small figure calibration

%!for other calibrations just copy move data and update format

%% Matlab start
close all
clear all
clc

%% Define global variables
global folder date cal rec test ctrl prop ...
    chkboard Cfit Dmap Pmap filesel

%% Additional Paths
dir= cd('..\'); % one folder up
addpath(genpath([cd '\Function_Libraries']))
cd(dir) %  back to original folder
clear dir

%% Locations
folder='D:\Desktop\OpenSourceCode';
%'C:\Users\Koen Muller\Desktop\Data_Processing';
%'D:\Desktop';
date='';
%'\2017_08_07_to_14';
%'\2019_01_24_to_30';
%'\2018_08_12_to_18';
%'\2018_01_03_to_10';
%'\calib_swaraj';
%'\2017_08_07_to_14';
%'\2018_01_03_to_10';
cal='\Calibration';
%'\Oceanium_Calibration';
%'\calib_imdata';
%'\Oceanium_Calibration'; 
rec='\Calibration_Route_003';
%'\2017_08_14_Time_08p30am_Calibration_010'; % Recording
%'\2019_01_25_time_15p30_Distortion_Correction20_foc_inf_ap_2p0_fov_x000y000_ex_5e4mus_fr_1Hz';
%'\2019_01_25_time_15p30_Calibration_Route64_foc_inf_ap_2p0_fov_x000y000_ex_5e4mus_fr_1Hz';
%'\2018_08_17_time_15p00_distortion_calibration19_fo_inf_ap_2p8_fov_x000y000_ex_10e4_fr_01p00hz';
%'\2017_08_14_Time_08p30am_Calibration_010'; % Recording
%'\2018_08_17_time_15p00_distortion_calibration19_fo_inf_ap_2p8_fov_x000y000_ex_10e4_fr_01p00hz';
%'\2018_01_10_14p00_distortion_correction09_fo_inf_ap_2_fov_x000y000_ex_10e4us_fr_01p00hz';
%'\2018_01_10_14p00_distortion_correction11_fo_inf_ap_2_fov_x000y000_ex_10e4us_fr_01p00hz';
%'\exp_20_frm';
%'\2017_08_14_Time_08p30am_Calibration_010'; % Recording
test='\TestImage';
%'\2017_08_13_Time_13p20pm_Extra_Light_All_Donut_Escape_ET_15e3_mus_FR_40fps'; % Recording
%'\2017_08_13_Time_13p20pm_Extra_Light_All_Donut_Escape_ET_15e3_mus_FR_40fps'; % Recording
%'\Oceanium_Calibration\2019_01_25_time_15p30_Calibration_Route01_foc_inf_ap_2p0_fov_x000y000_ex_5e4mus_fr_1Hz';
%'\Oceanium_Calibration\2018_08_17_time_15p00_calibration_route01_fo_inf_ap_2p8_fov_x000y000_ex_10e4_fr_01p00hz';
%'\Oceanium_Calibration\2018_01_10_14p00_distortion_correction11_fo_inf_ap_2_fov_x000y000_ex_10e4us_fr_01p00hz';
%'\test_exp';
%'\2017_08_10_Time_21p20pm_Donut_Escape_01_ET_15e3_mus_FR_40fps';
%'\2018_01_09_18p05_small_mill_fo_inf_ap_2_fov_x480y660_ex_20e3us_fr_41p56hz';

%% create path for processing function storage
if exist([folder date cal vsl 'mfunc'],'dir')~=7
    mkdir([folder date cal vsl 'mfunc'])
end
addpath([folder date cal vsl 'mfunc'])

%% Controls
if 0
    [prop,ctrl] = proc_set; % general processing settings
else
    load([folder date cal vsl 'prop.mat'])
    load([folder date cal vsl 'ctrl.mat'])
end

%% create checkerboard
if 0
    chkboard=def_chkboard;
else
    load([folder date cal vsl 'chkboard.mat'],'chkboard')
end

%% symbolic mapping and curves
if 0
    [ Cfit ] = def_cfit;
    [ Dmap ] = def_dmap;
    [ Pmap ] = def_pmap;
else
    load([folder date cal vsl 'Cfit.mat'])
    load([folder date cal vsl 'Dmap.mat'])
    load([folder date cal vsl 'Pmap.mat'])
end

%% Part I - Perform the image processing over input image data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run the image processing script to %%%
%%% inspect and process the images.    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get data from images
if 0
    improc_chkboard % check number with param chkb.
    %close all
end

%% Plot image processing
if 0 
    plot_imageprocessing % all
end

%% Part II - Perform the calibration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run the consecutive scripts to calibrate. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% select board to process for distortion and calibration
if 0
    filesel=sel_chkboard; % customize inside
else
    load([folder date cal vsl 'filesel.mat'])
end

%% distortion correction
if 0
    rec_sinview
    %close all
end

%% single view calibration
if 0
    cal_sinview
    %close all
end

%% multiple view calibration
if 0
    cal_mulview
    %close all
end

%% evaluate image warping
if 0
    def_cmap
end

%% statistics
if 0
    exp_stats
end

%% Part III - Plot results and use the calibration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot results corresponding to the calibration %%%
%%% , use the calibration by clicking the images. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot Calibration
if 0
    plot_calibrationimage(1,1,1) % input: file, camera, frame
end
if 0
    plot_checkerboards
end

%% quality assessment
if 0
    plot_errordistribution % not publish
end
if 0
    plot_dispvec % not publish 
end
if 0
    plot_spatialsampling % publish
end
if 0
    plot_spatialaccuracy % publish
end

%% perform simple example measurement
if 0
    click_raytracing(1) % input: frame
end

%% sandbox

