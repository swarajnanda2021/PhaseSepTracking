%% Info
%	
%   $Author: Koen Muller$
%
%   Functions to Patch background image
%
%-- Todo
%       All code have implemented variable slash for operating system but
%       need to be checked and where needed debugged

%% Matlab start
close all
clear 
clc

global folder date rec ctrl prop plotting

%% Additional Paths
dir= cd('..\'); % one folder up
addpath(genpath([cd '\Function_Libraries']))
cd(dir) %  back to original folder
clear dir

%% Locations
folder='D:\Desktop';%
%'C:\Users\Koen Muller\Desktop\Data_Processing\';
%'D:\Desktop';
%'G:\2016_Master_Thesis_Chlamydomonas_Swimming'; % Desktop Location
date='\2017_08_07_to_14';
%'\2018_08_12_to_18';
%'\2017_08_07_to_14';
%'\2018_01_03_to_10';
rec='\2017_08_10_Time_19p00pm_Milling_ET_40e3_mus_FR_20fps'; 
%'\2017_08_10_Time_20p40pm_Donut_Shape_Escape_ET_15e3_mus_FR_40fps';
%'\2018_01_09_18p05_small_mill_fo_inf_ap_2_fov_x480y660_ex_20e3us_fr_41p56hz';
%'\2019_01_29_time_21p50_Two_groups_interact_foc_inf_ap_2p0_fov_320y300_ex_40e3mus_fr_21p11Hz';
%'\2018_08_14_time_20p30_mill_split_foraging_fo_inf_ap_2p8_fov_x800y700_ex_40e3_fr_22p88hz';
%'\2018_08_16_time_14p20_foraging_group_fo_inf_ap_2p8_fov_x640y500_ex_20e3_fr_39p17hz';
% '\2018_08_16_time_22p55_milling_ball_fo_inf_ap_2p8_fov_x640y400_ex_20e3_fr_37p81hz';
%'\2017_08_13_Time_13p20pm_Extra_Light_All_Donut_Escape_ET_15e3_mus_FR_40fps'; % Recording
%'\2018_01_07_11p40_ball_milling_fo_inf_ap_2_fov_x480y400_ex_20e3us_fr_37p81hz';
%'\2017_08_10_Time_21p20pm_Donut_Escape_01_ET_15e3_mus_FR_40fps'; % Recording
%'\2018_01_09_18p05_small_mill_fo_inf_ap_2_fov_x480y660_ex_20e3us_fr_41p56hz';

%% Controls 
plotting='off'; % switch of all figures for cluster! ( 'off' | 'on')

[prop,ctrl] = proc_set; % processing setting function

%% Get data from images, or generator
if 0
    proc_background
end
