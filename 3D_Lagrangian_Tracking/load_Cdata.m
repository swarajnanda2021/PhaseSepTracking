function [Cpln,frm_new] = load_Cdata(frm_rem,frm_new)
%load_Cdata Load camera data from processed images with objects represented
%   as bounding ellipsoids over multiple scales.
%   
%   Input,
%       frm_rem Frames to remove from window data
%       frm_new Frames to append to window data
%   
%   Output,
%       Cpln Camera plane data
%   

%% Get global variables
global folder date rec 

%% load data
if exist([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Cpln.dat'],'file')~=0 % load
    
    % load previous image data
    Cpln=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Cpln.dat']);
    
    % trash processed frames previous data
    Cpln=Cpln(:,~ismember(Cpln(3,:),frm_rem));
    
    % remove processed frame if already passed this step
    frm_new=frm_new(:,~ismember(frm_new,Cpln(3,:)));
    
else % initiate
    
    % start image data
    Cpln=zeros(14,0); % [zero-index camera frame ellipse velocity intensity]
    
end

%% get new data there

% new camera data from 2D tracking
Cpln_new=fload([folder date rec vsl 'Preproc2DTracking' vsl 'Cpln.dat']); % better do memory map?

% append
Cpln=cat(2,Cpln,Cpln_new(:,ismember(Cpln_new(3,:),frm_new)));

end

