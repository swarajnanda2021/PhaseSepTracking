function postproc_vortextrajs



% Identify, join and smoothen tracks by fitting polynomial

% Track joining scheme:
% 3 and 4  frame tracks => extended by one timestep max, linear extrapolation
% 5 and 6  frame tracks => extended by two timesteps max, quadratic
% 7 and 8  frame tracks => extended by two timesteps max, quadratic
% 9 onward frame tracks => extended by three timesteps max

% JoinScheme = [  3 4 1 1;...
%                 5 8 2 2;... 
%                 9 0 3 2];           % it says from 3 to 4, extend by 1; from 5 to 8, extend by 2; from 9 onwards (0 being onwards), extend by 3 
            
            % JoinScheme's 4th column is the order of curve chosen
Lmin = 4;                           % Minimum track length that can be considered by this code 
Lmax = 50;                          % Maximum track length that can be considered by this code

% get globals
global folder date rec prop ctrl post traj

%% Load Track data
Pobj=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Pobj.dat']);


%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'VortexTrajectories'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'VortexTrajectories'])
end


%% remove nans, infs and tracks below Lmin
disp('post proc')
tic

%-- #) remove nans from preprocessing
N=~isnan(sum(Pobj,1)) & ~isinf(sum(Pobj,1));

%-- #) remove tracks smaller than Lmin timesteps
l=unique(Pobj(1,:));
cnt=histcounts(Pobj(1,:),[l-1/2,max(l)+1/2]);
%figure; histogram(cnt,[0.5:1:50.5])
% L=ismember(Pobj(1,:),l(cnt>=Lmin)); % at least a second/ characteristic timescale
L=ismember(Pobj(1,:),l(cnt>8 & cnt<10 )); % for testing

% valid data
val=N & L ;

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
%%
accum_disp=[];
for i=randperm(length(Posdata_ind_unq)) % 7,15,21 %15244
%     fprintf('%d ', i);
    i
    badtrack=0; % reset badtrack
    % find column index
    track_col_ind = find(Posdata(1,:)==Posdata_ind_unq(i));
    
    
    
    % Track length
    track_len = length(track_col_ind);
    % Store track data in dat format viz. rowwise is ones,time-step,x,y,z
    dat = [ones(size(Posdata(1,track_col_ind))); Posdata(2:5,track_col_ind)];
    
    
    
    % Begin RANSAC
    outlier = 1;     % Initiate outlier>inlier
    inlier = 0;
    dist_seg    = 1.8e-1; % maximum euclidean distance from fit to point 
    iter=1;
    itermax = 15;
    iterrestart = 1;
    nPosOutlier=0; % possible outliers in track
    multiplier = 0;
    vel_upper_outlier = 15;
    vel_lower_outlier = 1;
    while iter<itermax
%         iter
        if iterrestart==1 % Randomly sample track_len-1 (assuming one outlier) points and keep adding if they start meeting conditions
           temp_indx = randperm(track_len); 
           col_indx = temp_indx(1:track_len-nPosOutlier);
           dat2     = dat(:,col_indx);
        end
        
        if size(dat2,2)==3
%             ords = 1;
            ords = 2; % use quadratic always for smaller tracks (linear is problemmatic)
        elseif size(dat2,2)>3 && size(dat2,2)<6
            ords = 2;
        else
            ords = 2;
        end
        % Fit linear/quadratic/cubic model
        [ coef ,~] = polytraj( dat2, ords , 0 ); % estimate residual over all 
        
        % Calculate polynomial coefficient over subpoints
        [ tdat ] = polyeval( coef, dat, 0 );
        % Estimate euclidean residual over the entire track
        EucRes = vecnorm(abs(tdat - dat(3:5,:)));
        % Calculate velocity of each point using fit
        vel_dat = vecnorm(polyeval( coef, dat, 1 )).*prop.fps.*1e-3; % Calculate vel
        acc_dat = vecnorm(polyeval( coef, dat, 1 )); % Calculate double gradient
        
        
        % Distance segregation
        outlier_indx = find(EucRes>=dist_seg | abs(vel_dat)>=vel_upper_outlier | abs(vel_dat)<=vel_lower_outlier | isoutlier(acc_dat) ); % some thresholds for outlier
        inlier_indx = 1:size(dat,2);
        inlier_indx(find(ismember(inlier_indx,outlier_indx))) = [];
        inlier  = length(inlier_indx);
        outlier = length(outlier_indx);
        
        
        
        % Is outlier>inlier? 
        if outlier > floor((1+multiplier)*inlier) % If yes, restart iterations
            iterrestart = 1;
            
        else % Append inliers to dat2
            dat2 = dat(:,inlier_indx);
            iterrestart=0;
            
        end
        
        
        
        iter = iter+1; % CONTINUE ITERATIONS REGARDLESS
        
%         figure(1)
%         hold all
%         plot3(Posdata(3,track_col_ind),Posdata(4,track_col_ind),Posdata(5,track_col_ind),'k')
%         scatter3(Posdata(3,track_col_ind),Posdata(4,track_col_ind),Posdata(5,track_col_ind),[],vecnorm(polyeval( coef, dat, 1 )))
%         scatter3(tdat(1,:),tdat(2,:),tdat(3,:),[],vecnorm(polyeval( coef, dat, 1 )),'filled')
%         plot3(tdat(1,:),tdat(2,:),tdat(3,:),'b')
%         scatter3(tdat(1,inlier_indx),tdat(2,inlier_indx),tdat(3,inlier_indx))
%         axis equal
%         axis tight
%         drawnow 
%         pause
%         clf
%         
        if (iter==itermax) && (outlier>floor((1+multiplier)*inlier)) && (track_len>6) && (nPosOutlier<3) 
            iter=1;
            nPosOutlier=nPosOutlier+1;% increase number of outliers max to 3 outliers
%             disp('noutliers possible incremented')
        end
        
        if nPosOutlier>floor(0.5*track_len) && (ords==1) % if more than half the track are possibly outliers, just get rid of the track (can be problematic for small tracks). Also, if only a 1st order fit is possible, disregard it
           iter=itermax; 
           badtrack=1;
        end
        
        
        
        
        
    end
    
    
    
    % Store residual
    Res = [Res EucRes];
    
    % Store coefficients of the polynomial
    coef(1,:) = coef(1,:).*Posdata_ind_unq(i);
    Coef_mat = [Coef_mat coef];
    % Append to PosObj matrix with filtered tracks
    dat(1,:) = dat(1,:).*Posdata_ind_unq(i);
    PosObj = [PosObj dat(:,inlier_indx)];
    
    % Estimate 2 positions projected forward and backward in time and store
    % them in PosObj_predict
    nPredict=2;
%     t_ind = min(dat(2,:))-nPredict:max(dat(2,:))+nPredict;
    t_ind = min(dat(2,inlier_indx)):0.1:max(dat(2,inlier_indx));
    t_ind = [ones(size(t_ind)).*dat(1,1);t_ind];
    [ tdat ] = polyeval( coef, t_ind, 0 );
    
    
    
%     if outlier>floor((1+multiplier)*inlier)
%        disp(['Bad track with ' num2str(outlier) ' outliers and ' num2str(inlier) ' inliers.'])
%     else
%        disp(['Good track with ' num2str(outlier) ' outliers and ' num2str(inlier) ' inliers.']) 
%     end

    plotting=1;
    if plotting==1%(badtrack==1) && (floor((1+multiplier)*inlier)>outlier) && (plotting==1) %&& (outlier>0)
        
        figure(2)
%         subplot(1,2,1)
        view(2)
        title(['i/outlier/inlier/meanRes/maxRes/Vel:::: ' num2str(i) '/' num2str(outlier) '/' num2str(inlier) '/' num2str(mean(EucRes)) '/' num2str(max(EucRes)) '/' num2str(mean(vel_dat))])
        hold all
        plot3(Posdata(3,track_col_ind),Posdata(4,track_col_ind),Posdata(5,track_col_ind),'r--')
        scatter3(tdat(1,:),tdat(2,:),tdat(3,:),[],vecnorm(polyeval( coef, t_ind, 1 )).*prop.fps.*1e-3,'filled')
        scatter3(dat(3,outlier_indx),dat(4,outlier_indx),dat(5,outlier_indx),'r')
        scatter3(dat(3,inlier_indx),dat(4,inlier_indx),dat(5,inlier_indx),'b')
%         scatter3(tdat(1,inlier_indx),tdat(2,inlier_indx),tdat(3,inlier_indx),'k')
%         scatter3(tdat(1,outlier_indx),tdat(2,outlier_indx),tdat(3,outlier_indx),'rX')
%         plot3(tdat(1,inlier_indx),tdat(2,inlier_indx),tdat(3,inlier_indx),'k')
%         scatter3(Posdata(3,track_col_ind),Posdata(4,track_col_ind),Posdata(5,track_col_ind),[],vecnorm(polyeval( coef, dat, 1 )).*prop.fps.*1e-3,'.')
%         scatter3(tdat(1,inlier_indx),tdat(2,inlier_indx),tdat(3,inlier_indx),[],vecnorm(polyeval( coef, dat(:,inlier_indx), 1 )).*prop.fps.*1e-3,'filled')
%         scatter3(tdat(1,outlier_indx),tdat(2,outlier_indx),tdat(3,outlier_indx),[],vecnorm(polyeval( coef, dat(:,outlier_indx), 1 )).*prop.fps.*1e-3)
%         axis equal
%         axis tight
        colorbar
%         zlim([-20 100])
        caxis([3 8])
        colormap(jet)
        
        accum_disp = [accum_disp vecnorm([Posdata(3,track_col_ind(2:end)) - Posdata(3,track_col_ind(1:end-1));Posdata(4,track_col_ind(2:end)) - Posdata(4,track_col_ind(1:end-1));Posdata(5,track_col_ind(2:end)) - Posdata(5,track_col_ind(1:end-1))])./(Posdata(2,track_col_ind(2:end))-Posdata(2,track_col_ind(1:end-1)))];
        
%         subplot(1,2,2)
% %         hold all
%         histogram(accum_disp,0.3:0.1:5)
% %         axis equal
%         
%         drawnow 
        
        
        pause
        clf
    end
    
    
end

























end

