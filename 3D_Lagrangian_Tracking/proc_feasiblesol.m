function proc_feasiblesol(frm_rem,frm_proc,frm_new)
%proc_feasiblesol Process feasible solution which fits optimal to the image 
%   data over the processed data span.
%
%   Structure: Temporal loop with sub scripts.
%   
%   Input:
%       frm_rem     Frame removal
%       frm_proc    Processing frames
%       frm_new     New frames to process.
%
%   Load data:
%       Cpln    Indexed camera data
%       Tlink   Indexed temporal links
%       Mset    Indexed matched sets.
%
%   Output:
%       Fsol Feasible solution.
%
%   Working principle:
%       This code branches feasible solution starting from the image plane.
%      	We first backtrack to append new track segments and then forward,
%       tracks to continue track and connect them:
%           1) Seed new solutions and branching temporal links. 
%           2) Branch matched sets and filter for temporal consistency by
%              interpolation between data frames
%			3) Split track violating the reprojection error
%           4) Evaluate the new branched solutions (smoothing)
%           5) merge connecting tracks
%           6) Optimize seeded solution for unique trajectory extrapolation
%              with minimal overlap. This means that seeded trajectory is
%              considered a unique object trajectory and cannot splits in
%              to which would compromise writing the physical object next.
%
%   Note: There are possible improvement to be made in the runtime of
%       eva_Fsol. brn_Fsol and eva_Fsol lend itself for kalman filters.
%       Backtracking goes in stable memory usage, forward tracking can 
%       experience some memory spikes in high [optical] density tracking.
%   

%% Get globals
global folder date rec ctrl

%% Load data
Cpln=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Cpln.dat']);
Tlink=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Tlink.dat']);

if ctrl.nparts>1000
    Mset_old=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Mset.dat']);

    Mset = obj_disp2(Mset_old); % change vector entries
else
    Mset=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Mset.dat']);
end
    
%% initiate data
if exist([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Fsol.dat'],'file')~=0 
    
    % solution 
    Fsol=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Fsol.dat']);
    
    % remove processed frame if already passed this step
    frm_new=frm_new(:,~ismember(frm_new,Fsol(5,:)));
    
else
    
    % Feasible camera tracking
    Fsol=zeros(20,0); % [feas-index track-index ident-index camera-index frame-index conic vel Pos Vel Reproj]
    
end

%% Loop Frames to process

% Sweep backward over new frames and forward over processing frames
for n=[frm_new fliplr(frm_proc) frm_proc] 
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initiate and seed new tracks %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if 1
        
        Fsol=ini_Fsol(Cpln,Mset,Fsol,n);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Branch feasible tracks in cameras and from object space %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if 1
        
        Fsol=brn_Fsol(Cpln,Tlink,Fsol,n);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if 1
        
        Fsol=spl_Fsol(Fsol,n);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Evaluate feasible solution %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if 1
        
        Fsol=eva_Fsol(Fsol,n);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Merge feasible solution from seeding and backtracking %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if 1
        
        Fsol=mer_Fsol(Fsol,n);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Optimize feasible solutions %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if 1
        
        Fsol=lopt_Fsol(Fsol,n);
        
    end
    
end % n

%% Remove unneccessary ghost tracks by global optimum

% global optimization over frm_proc
Fsol=gopt_Fsol(Fsol,frm_proc);

%% Remove previous processing step

% remove frames
Fsol=Fsol(:,~ismember(Fsol(5,:),frm_rem));

%% save feasible tracking data

% feasible track solution
fsave([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Fsol.dat'],Fsol,'w'); % write

end

