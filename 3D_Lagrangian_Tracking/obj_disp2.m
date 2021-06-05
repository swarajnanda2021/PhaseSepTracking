function [Mset_lamboseen] = obj_disp2(Mset)

% 19/09/2020
% Currently supports only solid-body rotation.
% This script should get more complicated once irrotational motion needs to
% be specified

% 22/09/2020
% Script upgraded to a lamb-oseen vortex with 2-step axis estimation using
% prediction vectors from the outer rim of the vortex. Parameters need to
% be added to function line instead of in-code declaration.

%% Initiate data

global folder date rec ctrl

disp('Adjusting Mset')
% Make Mset unique
Mset_lamboseen = [];

for tstep= min(Mset(3,:)):max(Mset(3,:))
    
    
    idxmset = Mset(3,:)==tstep;
    mset = Mset(:,idxmset==1);


    % segregate particles near the centre
    exclusion_width = -15;
    mset_outer = (mset(5,:)<exclusion_width) +  (mset(5,:)>-exclusion_width) + (mset(6,:)<exclusion_width) + (mset(6,:)>-exclusion_width);

    mset=mset(:,mset_outer==1 | mset_outer==2);

    if 0 %plot check
        % plot triangulation
        figure(10); hold all;plot3(mset(5,:),mset(6,:),mset(7,:),'.'); axis equal; axis tight;view(2); drawnow
    end
   
    
    %% Estimate the direction of tangential motion
    
    % We estimate here the direction of tangential modtion of the
    % particles. We take the prediction vectors from the local 2D
    % correlation that is the output of opt_Cpln, and, included in Mset,
    % remove the points in the middle of the measurement domain, subsample
    % it, extract rotation and translation from a decomposition of the
    % displacement vector, remove outliers, redo the linear decomposition
    % again, yielding a decent estimate of the direction of the swirl, and,
    % the linear displacement.


    %% Step 1: Take sample of points for approximate axis estimation

    

    % plot(mset(5,:)<-exclusion_width)

    % take subsample of mset
    subsample = 2000;
    subsample_idx = randi([1 size(mset,2)],subsample,1);
    sub_mset = mset(:,subsample_idx);
    
    
    if 0 % plot check
        % plot triangulation
        figure(1); plot3(sub_mset(5,:),sub_mset(6,:),sub_mset(7,:),'.'); axis equal; view(2); drawnow

        % plot triangulation
        hold on
        % quiver3(sub_mset(5,:),sub_mset(6,:),sub_mset(7,:),sub_mset(8,:),sub_mset(9,:),sub_mset(10,:),0,'.')
        plot3(sub_mset(5,:)+sub_mset(8,:),sub_mset(6,:)+sub_mset(9,:),sub_mset(7,:)+sub_mset(10,:),'.'); axis equal; view(2); drawnow
        hold off
        axis tight

    end
    
    
    % perform decomposition into rotation, translation and scaling
    % components
    X_sub = [sub_mset(5,:) ; sub_mset(6,:) ; sub_mset(7,:)]'; % reconstructed 3D points, subsample
    Y_sub = [sub_mset(5,:)+sub_mset(8,:) ; sub_mset(6,:)+sub_mset(9,:) ; sub_mset(7,:)+sub_mset(10,:)]'; % estimated displacement in the next frame, subsample

    
    [d,Z,tr] = procrustes(Y_sub,X_sub); % Procustes analysis


    %% Step 2: Estimating the outliers and inliers based on distance threshold in the outer region of the axis and re-estimating a better axis
    % Store points and displacement vector in matrices
    mset = Mset(:,idxmset==1);
    X_all = [mset(5,:) ; mset(6,:) ; mset(7,:)]'; % reconstructed 3D points, full sample
    Y_all = [mset(5,:)+mset(8,:) ; mset(6,:)+mset(9,:) ; mset(7,:)+mset(10,:)]'; % estimated displacement in the next frame, full sample

    % Store linear transformation information    
    c = (ones([size(X_all,1) 3]).*tr.c(1,:)); % set translation vector
    T = tr.T;                                 % set rotation vector
    b = tr.b;                                 % set scale

    % Convert rotation matrix to axis-angle representation to extraction
    % principle rotation axis
    axang_rep = rotm2axang(tr.T);  % to save later or maybe never :P
    
    
    % Replace first three elements of axang_rep, i.e., the axis vector,
    % with the one retrieved from the volume carving section
    
     % Load visual hull and estimate axis
    load([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(tstep) '.mat'],'Vox_int','X_vox', 'Y_vox', 'Z_vox')
    % Reshape 3D matrix to 2D with transverse dimensions squashed    
    reshp_Vox_int = reshape(Vox_int,size(Vox_int,1)*size(Vox_int,2),size(Vox_int,3)); % Intensity
    reshp_X_Vox = reshape(X_vox,size(Vox_int,1)*size(Vox_int,2),size(Vox_int,3)); % X-coordinate
    reshp_Y_Vox = reshape(Y_vox,size(Vox_int,1)*size(Vox_int,2),size(Vox_int,3)); % Y-coordinate
    
    % Add X and Y coordinate matrices
    Centroid_cav = [ [sum(reshp_X_Vox.*reshp_Vox_int,1) ; sum(reshp_Y_Vox.*reshp_Vox_int,1)]./sum(reshp_Vox_int,1) ; squeeze(Z_vox(1,1,:))'];
    
    
    
    % Calculate unit vector from first and last points
    new_axis = abs(Centroid_cav(:,1) - Centroid_cav(:,end))./norm(Centroid_cav(:,1) - Centroid_cav(:,end));
    % Replace in axang_rep
    axang_rep(1) = new_axis(2);
    axang_rep(2) = new_axis(1);
    axang_rep(3) = new_axis(3);
    
    
    
    % Use axis representation to estimate shortest distance
    % For each point's 1X3 position <P_i> to origin <O> = [0 0 0] and axis
    % vector <U> [axang_rep(1) axang_rep(2) axang_rep(3)], the shortest
    % distance is given by:
    % d_i = abs(cross(<A_iO>,<U>))/abs(<U>);
    % where <> is the vector representation
    
    U = axang_rep(1:3); % storing the axis unit vector into U
    
    
    
    
    d = vecnorm(cross(X_all,repmat(U,size(X_all,1),1)),2,2)./vecnorm(U,2,2); % vector operation performing the distance calculation relative to this axis
    
    % Transform all 3D points to the new coordinate system using the
    % estimated linear transformation
    Z_all = b*X_all*T + c;

    % Estimate outliers in the outer region of the flow, 25+/-1 mm from the
    % current axis
 
    d_thresh_mid = 15;
    d_thresh_spread = 1;
    dist_thresh = 0.06; % outlier distance threshold
    testidxs = (d<(d_thresh_mid+d_thresh_spread)) & (d>(d_thresh_mid-d_thresh_spread)); % performing the radial distance thresholding
        
    % Store the outer-rim 3D particles into matrices
    X_test = X_all(testidxs==1,:);
    Y_test = Y_all(testidxs==1,:);
    Z_test = Z_all(testidxs==1,:);
    % Calculate difference between Koen's vectors and linear transformation
    % vectors at the outer-rim region of the present rotation axis
    dist = vecnorm(abs(Z_test(:,1:2)-Y_test(:,1:2)),2,2);
    
    % Store new candidates for rotation axis estimation in a different
    % matrix
    distthreshidx = dist<dist_thresh; % Performing the difference in prediction thresholding
    X_sub_2 = X_test(distthreshidx==1,:); % storing into matrixes
    Y_sub_2 = Y_test(distthreshidx==1,:);
    
    % Re-estimating the rotation matrix using this set of points
    [d_2,Z_sub_2,tr_2] = procrustes(Y_sub_2,X_sub_2); % Procustes analysis
    
    
    %% Estimating the irrotational vortex model
    
    % Estimate new transformation of Z_all using the latest global rotation
    % matrix
    Mset_new = [];
    c = (ones([size(X_all,1) 3]).*tr_2.c(1,:));
    Z_all = tr_2.b*X_all*tr_2.T + c;
    
    % Some parameters
    r_v = ctrl.rv*1e-3; % viscous core radius
    r_est = 25e-3; % outer-rim radius from global axis, where the lamb oseen will be fitted to
    f_aq = 5000; % camera acquisition frequency (can be changed later on, too lazy)
    delta_t = f_aq^-1;
    axang_rep = rotm2axang(tr.T);  % convert rotation matrix to axis angle, replace this value
    theta_est = axang_rep(4); % the value of angular displacement at r_est
    
    % vortex model is u_theta = (lambda_inf/(2 pi r))(1 - exp(-r^2/r_v^2),
    % and, u_theta = r*omega_theta; omega_theta is the angular velocity,
    % Further, as we kno    w omega_theta at r=25mm, we know u_theta at 25 mm
    % u_theta(25mm) = r_est*omega_theta_estimated = r_est *theta_est/delta_t;
    
    % this can be substituted to get some estimate lambda_inf for the
    % lamb-oseen vortex 
    % lambda_inf =  (r_est *theta_est/delta_t) * (2 pi r_Est) * (1 - exp(-r_est^2/r_v^2)^-1 ;
    lambda_inf =  (r_est *theta_est/delta_t) * (2*pi*r_est) * (1 - exp(-r_est^2/r_v^2))^-1 ;
    
    % make function
    theta_fin = @(rspace) ((lambda_inf./(2.*pi.*rspace)).*delta_t.*(1 - exp(-rspace.^2/r_v^2)))./rspace;
    
    % Using the function, calculate a rotation matrix for each partile
    rotmatrices_all = axang2rotm([repmat(axang_rep(1:3), size(X_all,1),1)  theta_fin(d.*1e-3)]);
    
    % transform X_all using d and theta_fin
    
    for iii=1:size(X_all,1)
    
        Z_lamboseen(iii,:) = b*X_all(iii,:)*squeeze(rotmatrices_all(:,:,iii)) + c(iii,:);
        
    end
   
    % Store data
    mset_new = Mset(:,idxmset==1);
    mset_new(8:10,:) = Z_lamboseen'-mset(5:7,:);
    % append
    Mset_lamboseen = [Mset_lamboseen mset_new];
    
    if 0
        
        % plot triangulation
        figure; plot3(mset_new(5,:),mset_new(6,:),mset_new(7,:),'.'); axis equal; view(2); drawnow
        
        % plot triangulation
        hold on
        quiver3(mset_new(5,:),mset_new(6,:),mset_new(7,:),mset_new(8,:),mset_new(9,:),mset_new(10,:),0,'.')
        plot3(mset_new(5,:)+mset_new(8,:),mset_new(6,:)+mset_new(9,:),mset_new(7,:)+mset_new(10,:),'.'); axis equal; view(2); drawnow
        hold off
        
        pause
        
    end


    if 0
        figure(1)
        clf
        histogram(vecnorm(abs(Y_all-Z_all),2,2))
        hold all
        histogram(vecnorm(abs(Y_all-Z_lamboseen),2,2))
        
        figure(2)
        clf
        hold all
        plot3(X_all(:,1),X_all(:,2),X_all(:,3),'.')
        plot3(Y_all(:,1),Y_all(:,2),Y_all(:,3),'.')
        quiver3(X_all(:,1),X_all(:,2),X_all(:,3),Y_all(:,1)-X_all(:,1),Y_all(:,2)-X_all(:,2),Y_all(:,3)-X_all(:,3),0,'r.')

        plot3(Z_lamboseen(:,1),Z_lamboseen(:,2),Z_lamboseen(:,3),'.')
        quiver3(X_all(:,1),X_all(:,2),X_all(:,3),Z_lamboseen(:,1)-X_all(:,1),Z_lamboseen(:,2)-X_all(:,2),Z_lamboseen(:,3)-X_all(:,3),0,'k.')


        quiver3(0,0,0,-axang_rep(1),-axang_rep(2),-axang_rep(3),100,'k','LineWidth',2)
        
        xlim([-20 20])
        ylim([-20 20])
        zlim([0 40])
        
        
        view(2)
        
        axis tight
        axis equal
        drawnow
        pause
    end
    
    
    
    Z_lamboseen = [];
    X_all = [];
    Y_all = [];
    Z_all = [];
    
end






end

