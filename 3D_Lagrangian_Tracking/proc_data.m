function proc_data
%proc_data Process data in marging window.
%
%   Consequtive data processing:
%       proc_newframes      Process new frames data
%       proc_feasiblesol    Process feasible solution
%       proc_physobject     Write PhysicalObject data.
%   

%% Get global variables
global folder date rec ctrl

%% initiate windowed time margin
if exist([folder date rec vsl 'Preproc3DTracking' vsl 'Pobj.dat'],'file')~=0
    
    % get end data
    Pobj=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Pobj.dat']);
    
    % loop parameters
    frm_begin=ctrl.tproc(1); % beginning frame
    frm_end=ctrl.tproc(1)+ctrl.tstep*floor(range(ctrl.tproc)/ctrl.tstep); % ending frame   
    frm_step=max(Pobj(2,:))+1 : ctrl.tstep : frm_end; % frame stepping
    
else
    
    % loop parameters
    frm_begin=ctrl.tproc(1); % beginning frame
    frm_end=ctrl.tproc(1)+ctrl.tstep*floor(range(ctrl.tproc)/ctrl.tstep); % ending frame   
    frm_step=frm_begin : ctrl.tstep : frm_end; % frame stepping
    
end

%% loop window over temporal data
for n=frm_step
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% define time marging window %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    frm_win= n : min( n-1+ctrl.twin , ctrl.tproc(2) );
    if n==frm_begin % first frame step
        frm_new= frm_win ; % frames new in queue
        frm_proc= frm_win( 1 : ctrl.tstep ) ; % frames to process
        frm_rem= [] ; % no physical frames to remove
    elseif n== frm_end % last frame step
        frm_new= [] ; % no new frames in queue
        frm_proc= frm_win ; % remaining frames to process
        frm_rem= [] ; % keep frames to inspect optimization in last window
    else % all intermediate steps otherwise
        frm_new= frm_win( ctrl.twin+1-ctrl.tstep : end ) ; % frames new in queue
        frm_proc= frm_win( 1 : ctrl.tstep ) ; % frames to process
        frm_rem= n-ctrl.tstep : n-1 ; % frames to remove from step
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% process new data frames %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if 1
        
        proc_newframes(frm_rem,frm_new)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Process feasible solution %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if 1
        
        proc_feasiblesol(frm_rem,frm_proc,frm_new) 
        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Write physical object data %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if 1
        
        proc_physobject(frm_proc)
        
    end
    
end

end