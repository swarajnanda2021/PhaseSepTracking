function [prop,ctrl] = proc_set
%proc_set Processing settings
%   
%   Input: user define processing settings here
%   
%	Output:
%       prop.mat Properties
%       ctrl.mat Control

%% get globals
global folder date cal rec test

%% properties 
%(probably a better way to get all info)
prop.ext='tiff';
[~,prop.res]=import_frames({folder date test},prop.ext,1,1);
prop.roi=[1 1 prop.res(1:2)]; %[1+480 1+660 480+prop.res(1) 660+prop.res(2)]; %[1 1 prop.res(1:2)];
[vid,prop.res]=import_frames({folder date cal rec},prop.ext,1,1);
dum=class(vid);
prop.bit=str2double(dum(5:end));
prop.fps=2;

prop.flens=24; % 7.5; % 24; % focal length lens [mm]
prop.nref=1.0003/1.363; % relative refractive index (salt)water/air [-]
prop.areachip=16.6*14; % area chip [mm^2]

%% controls
ctrl.fproc=1:5;%[1:5:20,20]; % frames to process
ctrl.roi=inf*[-1 -1 1 1];%[1 1 prop.res(1:2)];%+[750 600 -750 -600];

ctrl.fobj=strel('disk',4,0); % filter domain
ctrl.ford=[2 2 0]; % image regression

ctrl.cord=2; % n-order curvefit
ctrl.pdis=0.075; % distortion percentage estimate from qual. insp. image

ctrl.dmap='vector-polymap'; 
% 'interface-cor'
% 'vector-polymap'; 
% 'interface-mod';
% 'radial-polymap';
ctrl.mord=3;
% polymap rational order,
% 1 is a linear image mapping - terms in general "excluded"
% 2 is insufficient for vector map
% 3 is reasonable
% 4 is heavy

ctrl.nnod=[5 4]; % # checkerboard nodes in x y
ctrl.tsiz=0.3; % size square tiles checkerboard [m]

ctrl.optl=10^-6; % optimization tolerance, iter till perc reduce object func

%% save
save([folder date cal vsl 'prop.mat'],'prop')
save([folder date cal vsl 'ctrl.mat'],'ctrl')

end