function [Tlink,frm_new] = load_Tlink(frm_rem,frm_new)
%load_Tlink Load indexed temporal links to camera data
%   
%   Input,
%       frm_rem Frames to remove from window data
%       frm_new Frames to append to window data
%   
%   Output,
%       Tlink Indexed temporal links between ellipse identification

%% Get global variables
global folder date rec 

%% check for available data create tspan and indexing

if exist([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Tlink.dat'],'file')~=0 % load
    
    % load previous timelink data
    Tlink=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Tlink.dat']);
    
    % trash
    Tlink=Tlink(:,~ismember(Tlink(5,:),frm_rem));
%     Tlink=Tlink(:,~ismember(Tlink(1,:),...
%         Tlink(1,~ismember(Tlink(4,:),Cpln(1,:))...
%         | ~ismember(Tlink(6,:),Cpln(1,:))))); % remove relations
    
    % remove processed frame if already passed this step
    frm_new=frm_new(:,~ismember(frm_new,Tlink(5,:)));
%     frm_new=frm_new(:,~ismember(frm_new,Cpln(3,...
%         ismember(Cpln(1,:),Tlink([4 6],:))) ));
    
else % initiate
    
    % start time link data
    Tlink=zeros(7,0); % [track camera frame(n-1) i-camera frame(n) j-camera avg-disp-ellipse]
    
end

%% Get data from load

% new camera data from 2D tracking
Tlink_new=fload([folder date rec vsl 'Preproc2DTracking' vsl 'Tlink.dat']); % memory map?

% append
Tlink=cat(2,Tlink,Tlink_new(:,ismember(Tlink_new(5,:),frm_new))); % ismember(Tlink_new(3,:),frm_new) & 

end

