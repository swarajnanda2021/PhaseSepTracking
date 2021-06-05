function postproc_vortextrajs2

% Code that characterizes particle trajectories around the cavity and
% performes outlier detection (not removal). Outputs a new Pobj.dat file
% that contains the characterization and the statistical score of the
% track.

% get globals
global folder date rec prop ctrl post traj

%% Load data
% Track data
Pobj=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Pobj.dat']);
% Simplified Visual Hull data
load([folder date rec vsl 'Postproc3DReco' vsl 'ReducedCavityRepresentation' vsl 'Visual_hull_reduced.mat'],'Visual_hull_reduced');

%% create path for processing function storage
if exist([folder date rec vsl 'Vortex3DTracking' vsl 'VortexTrajectories'],'dir')~=7
    mkdir([folder date rec vsl 'Vortex3DTracking' vsl 'VortexTrajectories'])
end



%% remove nans, infs and tracks below Lmin
disp('post proc')

% Set minimum length of track you would post-process
Lmin = 15;

%-- #) remove nans from preprocessing
N=~isnan(sum(Pobj,1)) & ~isinf(sum(Pobj,1));

%-- #) remove tracks smaller than Lmin timesteps
l=unique(Pobj(1,:));
cnt=histcounts(Pobj(1,:),[l-1/2,max(l)+1/2]);
%figure; histogram(cnt,[0.5:1:50.5])
L=ismember(Pobj(1,:),l(cnt>=Lmin)); % at least a second/ characteristic timescale
% L=ismember(Pobj(1,:),l(cnt==Lmin)); % for testing

% valid data
val=N & L ;

%% Loop over tracks to extract radius of curvature, outliers and track averaged angular velocity

% write tracking data into Tdata
Tdata=Pobj( : ,  val ); % tracking

% Convert Tdata to PosData
Q=Tdata(3:12,:); % get quadric

% get position (shape and orientation are unimportant in track joining)
[ X, ~, ~ ] = qvec2eshp( Q ); % convert quadric to its position only

% Construct Posdata matrix
Posdata = [Tdata(1:2,:) ;X;Tdata(13:end,:)];


% extract indices of unique tracks
Posdata_ind_unq = unique(Posdata(1,:));
% Loop over each unique track
Res         = [];
Coef_mat    = [];
PosObj      = [];

% Begin track index counter
track_ind_cnt = 1;
% This following matrix is the final product of this code
Pobj_2 = [];

% Begin looping over all tracks
for i=1:length(Posdata_ind_unq)
    disp(['Evaluating track number ' num2str(i)])
    %% Extract Track Data
    % find column index
    track_col_ind = find(Posdata(1,:)==Posdata_ind_unq(i));
    % Store track data in dat format viz. rowwise is ones,time-step,x,y,z
    dat = [ones(size(Posdata(1,track_col_ind))); Posdata(2:5,track_col_ind)];
    
    %% Remove track constituents outside the bounds of the visual hull axis 
    dat2 = dat(:,dat(5,:)> -10 & dat(5,:) < 50);
    % Track length
    track_len = length(dat2(1,:));
    % The following bit of code should only be considered if the track
    % length remains above 5 after the segmentation
    if track_len >=5
        %% Fit polynomial to track
        % Fit linear/quadratic/cubic model
        ords=3; % set polynomial order
        [ coef ,~] = polytraj( dat2, ords , 0 ); % estimate residual over all 
        % Calculate polynomial coefficient over subpoints
        [ tdat ] = polyeval( coef, dat2, 0 );
        % Estimate euclidean residual over the entire track
        EucRes = vecnorm(abs(tdat - dat2(3:5,:)));
        % Calculate velocity of each point using fit
        vel_dat = vecnorm(polyeval( coef, dat2, 1 )).*prop.fps.*1e-3; % Calculate vel

        %% Estimate the distance of trajectory constituents from the cavity centroid at its respective timestep

        % Loop over track constituents
        for k=1:track_len
            % Retrieve tstep
            tstep = dat2(2,k);
            % Retrieve Centroid information of visual hull
%             x_VH        = Visual_hull_reduced(1,:,tstep);
%             y_VH        = Visual_hull_reduced(2,:,tstep);
            y_VH        = Visual_hull_reduced(1,:,tstep); % try switching y and x 
            x_VH        = Visual_hull_reduced(2,:,tstep);
            orient_VH   = Visual_hull_reduced(3,:,tstep);
            MajAx_VH    = Visual_hull_reduced(4,:,tstep);
            MinAx_VH    = Visual_hull_reduced(5,:,tstep);
            z_VH        = Visual_hull_reduced(6,:,tstep);
            % Retrieve z-position of the trajectory constituent at tstep
            z_Tr = dat2(5,k);
            % Estimate from linear extrapolation the position of x and y of the
            % centroid at the z-locations of the trajectory constituents
            x_Tr        = interp1(z_VH,x_VH,z_Tr);
            y_Tr        = interp1(z_VH,y_VH,z_Tr);
            orient_Tr   = interp1(z_VH,orient_VH,z_Tr);
            MajAx_Tr    = interp1(z_VH,MajAx_VH,z_Tr);
            MinAx_Tr    = interp1(z_VH,MinAx_VH,z_Tr);
            % Estimate distance of trajectory constituent from the cavity axis
            % estimated at the z-position of the particle
            r_Tr        = sqrt((x_Tr - dat2(3,k))^2 + (y_Tr - dat2(4,k))^2 );
            % Append to global track matrix
            Pobj_2 = [Pobj_2 [track_ind_cnt tstep dat2(3:5,k)' x_Tr y_Tr orient_Tr MajAx_Tr MinAx_Tr r_Tr]'];
            
            
        end
        
        % Update the counter
        track_ind_cnt = track_ind_cnt + 1;
    
    end

    %% Plot check
    if 0 %Plot check function
        
        figure(1)
        hold all
        plot3(dat(3,:),dat(4,:),dat(5,:),'r+')
        plot3(tdat(1,:),tdat(2,:),tdat(3,:),'k')
        axis equal
        view(2)
        
    end
    
end



%% Save Data

save([folder date rec vsl 'Vortex3DTracking' vsl 'VortexTrajectories' vsl 'Pobj_2_Lmin_' num2str(Lmin) '.mat'],'Pobj_2'); % save vecgrid






end

