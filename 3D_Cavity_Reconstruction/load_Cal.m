function [Dmap,Pmap,Dpar,Kmat,Rmat,tvec,camprop]=load_Cal
%load_cal Load calibration from either processed mechanical calibration 
%   data Davis, or checkerboard calibration method and define projective
%   matrix for each camera

%% Get global variables
global folder date cal

%% Load data
Dmap=importdata([folder date cal vsl 'Dmap.mat']);
Pmap=importdata([folder date cal vsl 'Pmap.mat']);
Dpar=importdata([folder date cal vsl 'Dpar.mat']);
Kmat=importdata([folder date cal vsl 'Kmat.mat']);
Rmat=importdata([folder date cal vsl 'Rmat.mat']);
tvec=importdata([folder date cal vsl 'tvec.mat']);
camprop=importdata([folder date cal vsl 'camprop.mat'],'camprop');

end