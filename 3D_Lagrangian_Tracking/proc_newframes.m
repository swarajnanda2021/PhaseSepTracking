function proc_newframes(frm_rem,frm_new)
%proc_newframes Process new frames is the overhead function to process new
%   frames.
%
%   Structure: Consequtive scripts that input the new frames window
%
%   Consequtive subfunctions:
%       get_Cdata Ellipsoid identifications from camera data
%       ind_Tlink Index temporal links camera data
%       def_Mset  Define matched index sets between camera data
%
%   Note,
%       This function acts as a forcing to the next step
%       

%% Get global variables
global folder date rec

%% Process new data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform image processing %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist([folder date rec vsl 'Preproc2DTracking' vsl 'Cpln.dat'],'file')~=0
    
    [Cpln,frm_cpln]=load_Cdata(frm_rem,frm_new); % in case already did 2D_Object_Tracking
    
else
    
    [Cpln,frm_cpln]=get_Cdata(frm_rem,frm_new);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Index track segments and define temporal links %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist([folder date rec vsl 'Preproc2DTracking' vsl 'Tlink.dat'],'file')~=0
    
    [Tlink,frm_tlink]=load_Tlink(frm_rem,frm_new); % in case already did 2D_Object_Tracking
    
else
    
    [Tlink,frm_tlink]=ind_Tlink(Cpln,frm_rem,frm_new);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find correspondence matching between cameras %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist([folder date rec vsl 'Preproc3DMatching' vsl 'Mset.dat'],'file')~=0
   
	[Mset,frm_mset]=load_Mset(frm_rem,frm_new); % this is an manual hack when want to speed up reprocessing Fsol
	
else
    
    [Mset,frm_mset]=def_Mset(Cpln,frm_rem,frm_new);
    
end

%% Save data

% camera data
fsave([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Cpln.dat'],Cpln,'w'); % write

% temporal links
fsave([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Tlink.dat'],Tlink,'w'); % write

% matched sets
fsave([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Mset.dat'],Mset,'w'); % write

%% Write and append camera data output

% camera plane
fsave([folder date rec vsl 'Preproc3DTracking' vsl 'Cpln.dat'],Cpln(:,ismember(Cpln(3,:),frm_cpln)),'a'); % append new or remaining

% camera tracking
fsave([folder date rec vsl 'Preproc3DTracking' vsl 'Tlink.dat'],Tlink(:,ismember(Tlink(5,:),frm_tlink)),'a'); % append new or remaining

% matched sets
fsave([folder date rec vsl 'Preproc3DTracking' vsl 'Mset.dat'],Mset(:,ismember(Mset(3,:),frm_mset)),'a'); % append new or remaining

end
