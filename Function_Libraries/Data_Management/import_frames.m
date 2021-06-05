function [vid,res] = import_frames( location, ext , frames , cams )
%import_frames This file imports different types of movie frame files
%   Input,
%       location or folder / file cell input
%       extension avi | im7 | tif
%       frames
%       cams
%
%   $Author: Koen Muller$

switch ext % choose extension
    case 'avi'
        
        [vid,res]=import_avi([cell2mat(location) ...
            location{end}],frames);
        
    case 'mp4' % color channels kept - hard-coded
        
        [vid,res]=import_mp4([cell2mat(location(1:end-1)) ...
            vsl location{end}],frames);
        
    case 'im7'
        
        [vid,res]=import_im7([cell2mat(location)...
            vsl],frames,cams);
        
    case {'tiff','tif'}
        
        [vid,res]=import_tiff([cell2mat(location) ...
		    vsl] , frames, cams);
        
    case {'mat'}
        
        [vid,res]=import_mat([cell2mat(location) ...
		    vsl] , frames , cams );
        
end

end

%-- extension specific function file
function [vid,res]=import_avi(file,frames)
% resolution
dum=VideoReader([file '.avi']); % first frame to initiate data management
res=[dum.Width dum.Height dum.Duration*dum.FrameRate 1]; % denote there is only 1 camera in avi

% Read files
vid=zeros([res(1:2),length(frames)],...
    'like',readFrame(dum)); % initiate video data

for n=1:length(frames) % loop selected frames
    dum.CurrentTime=(frames(n)-1)/dum.FrameRate;
    vid(:,:,n)=flipud(readFrame(dum))'; % column dimensional array convention
end % n

end

function [vid,res]=import_mp4(file,frames) % color channels kept - hard-coded
% resolution
dum=VideoReader([file '.mp4']); % first frame to initiate data management
res=[dum.Width dum.Height dum.Duration*dum.FrameRate 1]; % denote there is only 1 camera in avi

% Read files
vid=zeros([res(1:2),3,length(frames)],...
    'like',readFrame(dum)); % initiate video data

for n=1:length(frames) % loop selected frames
    dum.CurrentTime=(frames(n)-1)/dum.FrameRate;
    vid(:,:,:,n)=permute(flipud(readFrame(dum)),[2 1 3]); % column dimensional array convention
end % n

end

function [vid,res]=import_im7(folder,frames,cams)
% content folder (wildcard)
cont=dir([folder '*.im7']); % directory contents

% resolution
dum=readimx([folder,cont(1).name]); % first frame to initiate data management
res=[size(cell2mat(dum.Frames{1}.Components{1}.Planes)),size(cont,1),size(dum.Frames,1)];

% Read im7 files
vid=zeros([res(1:2),length(frames),length(cams)],...
    'like',cell2mat(dum.Frames{1}.Components{1}.Planes)); % initiate video data

% loop
for n=1:length(frames) % loop selected frames
    dum=readimx([folder,cont(frames(n)).name]); % current frame
    for m=1:length(cams) % loop selected camera's
        vid(:,:,n,m)=fliplr(cell2mat(dum.Frames{cams(m)}.Components{1}.Planes)); % write
    end % m fliplr(if need to put im top true
end % n

end

function [vid,res]=import_tiff(folder,frames,cams)

% camera folders (wildcard)
camfol=dir([folder 'Camera_*.']); % assume impose folder structure

% dummy camera
frmfile=dir([folder camfol(1).name vsl '*.tiff']); % assume all folders same frame content

% resolution
dum = imread([folder camfol(1).name vsl frmfile(1).name]); % first frame to initiate data management
res=[size(dum,1),size(dum,2),length(frmfile),length(camfol)];

% Read tif/tiff files
vid=zeros([res(1:2),length(frames),length(cams)],...
    'like',dum); % initiate video data

% loop
for n=1:length(frames) % loop selected frames
    for m=1:length(cams) % loop selected cameras
        vid(:,:,n,m)=imread([folder camfol(cams(m)).name vsl frmfile(frames(n)).name ]); % write
    end % m
end % n

end

function [vid,res]=import_mat(folder,frames,cams)
% content folder (wildcard)
cont=dir([folder '*.mat']); % directory contents

% resolution
dum=load([folder,cont(1).name]); % first frame to initiate data management
dum=dum.Imgs;
res=[size(dum,2),size(dum,1),size(cont,1),size(dum,3)];

% Read mat files
vid=zeros([res(1:2),length(frames),length(cams)],...
    'like',dum); % initiate video data

% loop
for n=1:length(frames) % loop selected frames
    dum=load([folder,cont(frames(n)).name]); % current frame
    dum=dum.Imgs;
    for m=1:length(cams) % loop selected camera's
        vid(:,:,n,m)=fliplr(dum(:,:,cams(m))'); % write
    end % m fliplr(if need to put im top true
end % n

end