function postproc_CavityCondAvg




%% get globals
global folder date rec prop mset post%plotting

%% Load data
% load general data
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
% Time=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Time.dat']);
Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
load([folder date rec vsl 'Postproc3DTracking' vsl 'TrackQualityMetric.mat'],'t_spurious');

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'CavityReference'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'CavityReference'])
end



Rdata = zeros(1,0);
Uthdata = zeros(1,0);
    

% loop frames
for n=post.tproc(1):post.tproc(2)%mset.tspan(1):mset.tspan(2)
    disp(['Starting tstep ' num2str(n)])
    
    
    %% For visual hull processing
    % Load visual hull and estimate axis
    load([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(n) '.mat'],'Centroid_cav')
    % Reshape 3D matrix to 2D with transverse dimensions squashed    

        
    
    
    %% For velocity processing
    % Time set
    N=(ismember(Index(2,:),n));  
    outl = ~(ismember(Index(1,:),unique(t_spurious)));
    % Gather position and velocities
    pos = Position(:,N&outl);
    vel = Velocity(:,N&outl);
    
    % determine radial position
    centre_x = interp1(Centroid_cav(3,:),Centroid_cav(1,:),pos(3,:));
    centre_y = interp1(Centroid_cav(3,:),Centroid_cav(2,:),pos(3,:));
    
    r_pos = ((pos(1,:)-centre_x).^2 + (pos(2,:)-centre_y).^2).^0.5;
    theta_pos = atan2((pos(1,:)-centre_y),(pos(2,:)-centre_x));
    u_theta_pos = ((vel(1,:).^2) + (vel(2,:).^2)).^0.5;
    
    %% Gather data
    Rdata = cat(2,Rdata,r_pos);
    Uthdata = cat(2,Uthdata,u_theta_pos);
    
    
    %% save data at end of each timestep
    if 0
        disp(['Storing tstep ' num2str(n) ' completed!'])
        save([folder date rec vsl 'Postproc3DTracking' vsl 'CavityReference' vsl 'Tstep_' num2str(n) '.mat'],'pos_new','vec_accum_mat')
    end
    
    
end


%%


histogram2(Rdata,Uthdata./1000,[50 50],'Normalization','probability','FaceColor','flat');
% ylim([0 2])
% xlim([0 25])

view(2)




end

