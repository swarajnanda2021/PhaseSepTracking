%% Code for inverting the distortion mapping estimated from Koen's calibration code

close all
clear all
clc

global folder date rec cal ctrl prop plotting ...
    Kmat Rmat tvec Pmat Dpar Dmap Pmap camprop

%% Additional Paths
direc= cd('..\'); % one folder up
addpath(genpath([cd '\Function_Libraries']))
cd(direc) %  back to original folder
clear direc % remove variable

%% Locations
folder='C:\Work\PostGraduation\Scripts\ParticleTracking\Koen_testdata'; % Desktop Location
date='\Test_Swaraj'; % Date Folder
rec='\Experiment_Raw_23'; % Recording
cal='\Calibration_200514_141513'; % Calibration files

%% Load files
[Dmap,Pmap,Dpar,Kmat,Rmat,tvec,camprop]=load_Cal;

% choose frame
frame = 1;

% load distortion matrix
distortn = load([folder date cal vsl 'wrp_',num2str(frame),'.mat'],'wrp');

%% Calculate Dmap input output matrix

% make the input arrays for Dmap
imgsize = size(distortn.wrp.x);

pix_x = [1  : imgsize(1)];
pix_y = [1  : imgsize(2)];
[mgrid_pixx, mgrid_pixy] = meshgrid(pix_x,pix_y);

input_array_pixx = reshape(mgrid_pixx,1,length(pix_x)*length(pix_y));
input_array_pixy = reshape(mgrid_pixy,1,length(pix_x)*length(pix_y));

input_array_pixxy = [input_array_pixx; input_array_pixy];

% calculate the output arrays for Dmap
output_array_pixxy  = Dmap(Dpar{1,frame}(1),Dpar{1,frame}(2),Dpar{1,frame}(3),Dpar{1,frame}(4),Dpar{1,frame}(5),Dpar{1,frame}(6),Dpar{1,frame}(7),Dpar{1,frame}(8),Dpar{1,frame}(9),Dpar{1,frame}(10),Dpar{1,frame}(11),Dpar{1,frame}(12),Dpar{1,frame}(13),Dpar{1,frame}(14),Dpar{1,frame}(15),Dpar{1,frame}(16),Dpar{1,frame}(17),Dpar{1,frame}(18),Dpar{1,frame}(19),Dpar{1,frame}(20),input_array_pixx,input_array_pixy);



tform = fitgeotrans(output_array_pixxy',input_array_pixxy','polynomial',3);


%% Estimate mapping for converting output back to input
% 
% 
% % create minimization function
% 
% 
% x0 = rand(20,1);
% options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
% options.MaxFunctionEvaluations = 20000;
%     
% out = fminunc(@(param) inv_distortn(param,output_array_pixxy,input_array_pixxy),x0,options);
%     
% 
% %% Test the minimized functional
% params = out;
% % 
% dewarped = output_array_pixxy;
% 
% warped = input_array_pixxy;
% 
% 
% % Take cubic polynomial
% x = dewarped(1,:)';
% y = dewarped(2,:)';
% 
% 
% 
% x = dewarped(1,:)';
% y = dewarped(2,:)';
% 
% 
% 
% P = [ones(size(x)) x y x.*y x.^2 y.^2 ...
%          (x.^2).*y x.^3 (y.^2).*x ...
%         y.^3];
% 
% 
%     
% warped_est = [P*params(1:10,1) P*params(11:20,1)]' ;
% 
% 
% 
% residual = max(vecnorm((warped_est - warped)'));
% 



