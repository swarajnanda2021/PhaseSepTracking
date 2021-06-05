function [Mset,frm_new] = load_Mset(frm_rem,frm_new)
%load_Mset Load indexed correspondence matches
%   
%   Input,
%       frm_rem Frames to remove from window data
%       frm_new Frames to append to window data
%   
%   Output,
%       Mset Matched correspondence set
%
%	Note: Work similar to load_Cdata/Tlink only need
%		    to Hack: 
%		  Create a folder Preproc3DMatching and put 
%		    [opt_]Mset.dat/bin there to have this 
%		    function loading the data and not needing
%		    to process it in def_Mset.m.

%% Get global variables
global folder date rec 

%% check for available data create tspan and indexing

if exist([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Mset.dat'],'file')~=0 % load
    
    % load previous timelink data
    Mset=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Mset.dat']);
    
    % trash
    Mset=Mset(:,~ismember(Mset(3,:),frm_rem));
%     Tlink=Tlink(:,~ismember(Tlink(1,:),...
%         Tlink(1,~ismember(Tlink(4,:),Cpln(1,:))...
%         | ~ismember(Tlink(6,:),Cpln(1,:))))); % remove relations
    
    % remove processed frame if already passed this step
    frm_new=frm_new(:,~ismember(frm_new,Mset(3,:)));
%     frm_new=frm_new(:,~ismember(frm_new,Cpln(3,...
%         ismember(Cpln(1,:),Tlink([4 6],:))) ));
    
else % initiate
    
    % start time link data
    Mset=zeros(11,0); % [track camera frame(n-1) i-camera frame(n) j-camera avg-disp-ellipse]
    
end

%% Get data from load

% new camera data from 3D matching (manual created folder see above)
Mset_new=fload([folder date rec vsl 'Preproc3DMatching' vsl 'Mset.dat']); % memory map?

% append
Mset=cat(2,Mset,Mset_new(:,ismember(Mset_new(3,:),frm_new))); % ismember(Tlink_new(3,:),frm_new) & 

end

